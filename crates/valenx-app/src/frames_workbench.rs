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
use nalgebra::Vector3;

use valenx_frames::{cross_section_polygon, Profile};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

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

                    // Live cross-section outline (face-on).
                    if let Some(pts) = preview_polygon(s) {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Cross-section preview").strong());
                        draw_section_preview(ui, &pts);
                    }
                });
        });
}

/// Build the active [`Profile`] from the form and return its closed
/// cross-section outline as an XY polyline (z = 0), best-effort `None`
/// for an invalid spec. Drives the live face-on section preview and
/// re-validates the active family's dimensions like [`run_frames`].
fn preview_polygon(s: &FramesWorkbenchState) -> Option<Vec<Vector3<f64>>> {
    // Validate only the dimensions the active family uses — the same
    // preconditions as `run_frames`.
    let dims: &[f64] = match s.kind {
        ProfileKind::IBeam | ProfileKind::TBeam => &[s.h, s.b, s.tw, s.tf],
        ProfileKind::CChannel => &[s.h, s.b, s.tw],
        ProfileKind::LAngle | ProfileKind::RhsRect => &[s.h, s.b, s.t],
        ProfileKind::ChsRound => &[s.d, s.t],
    };
    if !dims.iter().all(|v| v.is_finite() && *v > 0.0) {
        return None;
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

    // `cross_section_polygon` returns a non-closed CCW outline; lift it
    // into the z = 0 plane and push the first vertex to close the loop.
    let outline = cross_section_polygon(profile);
    if outline.len() < 2 {
        return None;
    }
    let mut pts: Vec<Vector3<f64>> = outline
        .iter()
        .map(|p| Vector3::new(p[0], p[1], 0.0))
        .collect();
    let first = pts[0];
    pts.push(first);
    Some(pts)
}

/// Draw a closed XY polyline as a face-on (front-view, looking down −Z)
/// wireframe in a fixed-height canvas, the camera auto-framed to the
/// section's bounds.
fn draw_section_preview(ui: &mut egui::Ui, pts: &[Vector3<f64>]) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 200.0),
        egui::Sense::hover(),
    );
    let rect = response.rect;

    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            let v = p[k] as f32;
            min[k] = min[k].min(v);
            max[k] = max[k].max(v);
        }
    }

    let mut cam = OrbitCamera::default();
    cam.set_view(ViewDirection::Front);
    cam.frame_bounds(min, max);

    let (w, h) = (rect.width(), rect.height());
    let stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255));
    for pair in pts.windows(2) {
        let a = project_point(
            &cam,
            w,
            h,
            [pair[0].x as f32, pair[0].y as f32, pair[0].z as f32],
        );
        let b = project_point(
            &cam,
            w,
            h,
            [pair[1].x as f32, pair[1].y as f32, pair[1].z as f32],
        );
        if let (Some(a), Some(b)) = (a, b) {
            painter.line_segment(
                [
                    egui::pos2(rect.min.x + a.x, rect.min.y + a.y),
                    egui::pos2(rect.min.x + b.x, rect.min.y + b.y),
                ],
                stroke,
            );
        }
    }
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
        ProfileKind::CChannel => &[(s.h, "height"), (s.b, "width"), (s.tw, "web thickness")],
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

    #[test]
    fn preview_polygon_default_ibeam_is_closed_and_spans_h_by_b() {
        let s = FramesWorkbenchState::default();
        let pts = preview_polygon(&s).expect("default IPE 200 yields a section outline");
        // 12 I-beam vertices + the duplicated first vertex closing the loop.
        assert_eq!(pts.len(), 13);
        assert_eq!(pts.first(), pts.last());
        // Planar XY profile (z = 0) spanning b wide × h tall.
        assert!(pts.iter().all(|p| p.z == 0.0));
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for p in &pts {
            xmin = xmin.min(p.x);
            xmax = xmax.max(p.x);
            ymin = ymin.min(p.y);
            ymax = ymax.max(p.y);
        }
        assert!((xmax - xmin - s.b).abs() < 1e-9); // width  ≈ flange width b
        assert!((ymax - ymin - s.h).abs() < 1e-9); // height ≈ overall height h
    }

    #[test]
    fn preview_polygon_vertex_count_per_family() {
        // (family, expected closed-loop length = backend polygon verts + 1).
        for (kind, n) in [
            (ProfileKind::IBeam, 13),
            (ProfileKind::CChannel, 9),
            (ProfileKind::LAngle, 7),
            (ProfileKind::RhsRect, 5),
            (ProfileKind::ChsRound, 25),
            (ProfileKind::TBeam, 9),
        ] {
            let s = FramesWorkbenchState {
                kind,
                ..Default::default()
            };
            let pts = preview_polygon(&s)
                .unwrap_or_else(|| panic!("{kind:?} with default dims should preview"));
            assert_eq!(pts.len(), n, "{kind:?} closed-loop vertex count");
            assert_eq!(pts.first(), pts.last(), "{kind:?} loop is closed");
            assert!(pts.iter().all(|p| p.z == 0.0), "{kind:?} is planar (z = 0)");
        }
    }

    #[test]
    fn preview_polygon_none_for_invalid_active_dim() {
        // Zero height on the default I-beam → no preview (an inactive
        // family's dim being zero would still preview).
        let s = FramesWorkbenchState {
            h: 0.0,
            ..Default::default()
        };
        assert!(preview_polygon(&s).is_none());
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
