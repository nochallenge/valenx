//! The right-side **Frames Workbench** panel — structural cross-section
//! properties over `valenx-frames`.
//!
//! Mirrors the springs / gears / … / fasteners workbenches: a resizable
//! [`egui::SidePanel`] gated on [`crate::ValenxApp::show_frames_workbench`],
//! toggled from the View menu. The form picks a profile family (I-beam,
//! C-channel, L-angle, rectangular HSS, round CHS, or T-beam) and its
//! dimensions; the "Compute" button reports the cross-sectional area and
//! perimeter of the idealised (sharp-cornered) section, as a monospace
//! readout.

use eframe::egui;

use valenx_frames::Profile;

use crate::ValenxApp;

/// The profile family selected in the form. Mirrors the variants of
/// [`valenx_frames::Profile`] without the per-variant fields, which live
/// flat on the workbench state.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum ProfileKind {
    #[default]
    IBeam,
    CChannel,
    LAngle,
    RhsRect,
    ChsRound,
    TBeam,
}

/// Persistent form + result state for the Frames Workbench. Dimensions are
/// kept flat and reused across the profile families that share them.
pub struct FramesWorkbenchState {
    /// Selected profile family.
    kind: ProfileKind,
    /// Overall height/depth `h` (mm).
    h: f64,
    /// Overall width / flange width `b` (mm).
    b: f64,
    /// Web thickness `t_w` (mm) — I-beam / C-channel / T-beam.
    tw: f64,
    /// Flange thickness `t_f` (mm) — I-beam / T-beam.
    tf: f64,
    /// Wall / leg thickness `t` (mm) — L-angle / rect HSS / round CHS.
    t: f64,
    /// Outside diameter `d` (mm) — round CHS.
    d: f64,
    /// Formatted section readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for FramesWorkbenchState {
    fn default() -> Self {
        // An IPE 200 wide-flange (valenx_frames::Profile::default_ipe200),
        // plus sensible defaults for the other families' extra dimensions.
        Self {
            kind: ProfileKind::IBeam,
            h: 200.0,
            b: 100.0,
            tw: 5.6,
            tf: 8.5,
            t: 8.0,
            d: 100.0,
            result: String::new(),
            error: None,
        }
    }
}

/// One labelled `DragValue` row for a dimension.
fn dim_row(ui: &mut egui::Ui, label: &str, value: &mut f64) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(egui::DragValue::new(value).speed(0.5));
    });
}

/// Draw the Frames Workbench right-side panel. A no-op when the
/// `show_frames_workbench` toggle is off.
pub fn draw_frames_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_frames_workbench {
        return;
    }

    egui::SidePanel::right("valenx_frames_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            ui.heading("Frames");
            ui.label(
                egui::RichText::new("structural cross-section properties · valenx-frames")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.frames;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Profile").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.kind, ProfileKind::IBeam, "I-beam");
                        ui.radio_value(&mut s.kind, ProfileKind::CChannel, "C-channel");
                        ui.radio_value(&mut s.kind, ProfileKind::LAngle, "L-angle");
                        ui.radio_value(&mut s.kind, ProfileKind::RhsRect, "Rect HSS");
                        ui.radio_value(&mut s.kind, ProfileKind::ChsRound, "Round CHS");
                        ui.radio_value(&mut s.kind, ProfileKind::TBeam, "T-beam");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Dimensions (mm)").strong());
                    match s.kind {
                        ProfileKind::IBeam | ProfileKind::TBeam => {
                            dim_row(ui, "height h", &mut s.h);
                            dim_row(ui, "width b", &mut s.b);
                            dim_row(ui, "web t_w", &mut s.tw);
                            dim_row(ui, "flange t_f", &mut s.tf);
                        }
                        ProfileKind::CChannel => {
                            dim_row(ui, "height h", &mut s.h);
                            dim_row(ui, "width b", &mut s.b);
                            dim_row(ui, "web t_w", &mut s.tw);
                        }
                        ProfileKind::LAngle | ProfileKind::RhsRect => {
                            dim_row(ui, "height h", &mut s.h);
                            dim_row(ui, "width b", &mut s.b);
                            dim_row(ui, "thickness t", &mut s.t);
                        }
                        ProfileKind::ChsRound => {
                            dim_row(ui, "diameter d", &mut s.d);
                            dim_row(ui, "wall t", &mut s.t);
                        }
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_frames(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Section").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });
}

/// Build the selected [`Profile`] variant from the form, compute the
/// cross-sectional area + perimeter, and format the readout. Extracted
/// from the draw closure so it is unit-testable.
fn run_frames(s: &mut FramesWorkbenchState) {
    s.error = None;

    // Validate only the dimensions the selected family actually uses.
    let dims: &[(f64, &str)] = match s.kind {
        ProfileKind::IBeam | ProfileKind::TBeam => &[
            (s.h, "height"),
            (s.b, "width"),
            (s.tw, "web thickness"),
            (s.tf, "flange thickness"),
        ],
        ProfileKind::CChannel => {
            &[(s.h, "height"), (s.b, "width"), (s.tw, "web thickness")]
        }
        ProfileKind::LAngle | ProfileKind::RhsRect => {
            &[(s.h, "height"), (s.b, "width"), (s.t, "thickness")]
        }
        ProfileKind::ChsRound => &[(s.d, "diameter"), (s.t, "wall thickness")],
    };
    for (v, name) in dims {
        if !(v.is_finite() && *v > 0.0) {
            s.error = Some(format!("{name} must be positive"));
            return;
        }
    }

    let profile = match s.kind {
        ProfileKind::IBeam => Profile::IBeam {
            h: s.h,
            b: s.b,
            tw: s.tw,
            tf: s.tf,
        },
        ProfileKind::CChannel => Profile::CChannel {
            h: s.h,
            b: s.b,
            tw: s.tw,
        },
        ProfileKind::LAngle => Profile::LAngle {
            h: s.h,
            b: s.b,
            t: s.t,
        },
        ProfileKind::RhsRect => Profile::RhsRect {
            h: s.h,
            b: s.b,
            t: s.t,
        },
        ProfileKind::ChsRound => Profile::ChsRound { d: s.d, t: s.t },
        ProfileKind::TBeam => Profile::TBeam {
            h: s.h,
            b: s.b,
            tw: s.tw,
            tf: s.tf,
        },
    };

    s.result = format!(
        "profile            : {}\n\n\
         cross-section area : {:.2} mm\u{00B2}\n\
         perimeter          : {:.2} mm",
        profile.label(),
        profile.cross_section_area_mm2(),
        profile.cross_section_perimeter_mm(),
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = FramesWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_default_ibeam() {
        let mut s = FramesWorkbenchState::default();
        run_frames(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("cross-section area"));
        assert!(s.result.contains("perimeter"));
        assert!(s.result.contains("I-beam"));
        // IPE 200 (h=200, b=100, tw=5.6, tf=8.5): area = 2·100·8.5 + 5.6·183
        // = 2724.8 mm²; perimeter = 2·(200+200−5.6) = 788.8 mm. Recompute via
        // the backend to confirm the default.
        let p = Profile::IBeam {
            h: 200.0,
            b: 100.0,
            tw: 5.6,
            tf: 8.5,
        };
        assert!((p.cross_section_area_mm2() - 2724.8).abs() < 1e-9);
        assert!((p.cross_section_perimeter_mm() - 788.8).abs() < 1e-9);
        assert!(s.result.contains("2724.80"));
        assert!(s.result.contains("788.80"));
    }

    #[test]
    fn compute_round_chs_uses_diameter_and_wall() {
        let mut s = FramesWorkbenchState {
            kind: ProfileKind::ChsRound,
            d: 100.0,
            t: 4.0,
            ..Default::default()
        };
        run_frames(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("Round CHS"));
        // perimeter = π·d = π·100 ≈ 314.16 mm.
        assert!(s.result.contains("314.16"));
    }

    #[test]
    fn compute_rejects_nonpositive_active_dim() {
        // Zero height on the default I-beam → error (an inactive family's
        // dim being zero would NOT error).
        let mut s = FramesWorkbenchState {
            h: 0.0,
            ..Default::default()
        };
        run_frames(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
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
            draw_frames_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_frames_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_frames_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_frames_workbench = true;
        run_frames(&mut app.frames);
        app.frames.error = Some("invalid profile".to_string());
        draw_workbench(&mut app);
    }
}
