//! The right-side **Frames Workbench** panel — structural cross-section
//! properties over `valenx-frames`.
//!
//! Mirrors the springs / gears / … / fasteners workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_frames_workbench`,
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

impl FramesWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`.
    /// Covers every dimension caption across all profile families, so an agent
    /// can set a dimension regardless of the selected profile. `thickness t` and
    /// `wall t` both address the single `t` field (the caption differs by
    /// profile); `diameter d` is the round-CHS outside diameter.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "height h",
            "width b",
            "web t_w",
            "flange t_f",
            "thickness t",
            "wall t",
            "diameter d",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / non-positive value returns `Err(String)` — never a panic.
    /// Every cross-section dimension must be a finite `> 0` (mirrors the
    /// workbench's own `run_frames` validity check). The thickness field accepts
    /// both its `thickness t` (L-angle / rect HSS) and `wall t` (round CHS)
    /// captions.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let positive = |v: f64, what: &str| -> Result<f64, String> {
            if v.is_finite() && v > 0.0 {
                Ok(v)
            } else {
                Err(format!("{what} must be > 0, got {v}"))
            }
        };
        match name {
            "height h" => self.h = positive(value.as_f64()?, "height h")?,
            "width b" => self.b = positive(value.as_f64()?, "width b")?,
            "web t_w" => self.tw = positive(value.as_f64()?, "web t_w")?,
            "flange t_f" => self.tf = positive(value.as_f64()?, "flange t_f")?,
            "thickness t" | "wall t" => self.t = positive(value.as_f64()?, "thickness t")?,
            "diameter d" => self.d = positive(value.as_f64()?, "diameter d")?,
            other => return Err(format!("unknown Frames control: {other:?}")),
        }
        Ok(())
    }
}

/// One labelled `DragValue` row for a dimension.
fn dim_row(ui: &mut egui::Ui, label: &str, value: &mut f64) {
    ui.horizontal(|ui| {
        // Associate the `DragValue` with its caption via `labelled_by`, so the
        // spin button carries the caption as its accessibility / UI-Automation
        // Name (egui clears a DragValue's own name, so without this it is
        // anonymous to a screen reader / AI driver).
        let cap = ui.label(label);
        ui.add(egui::DragValue::new(value).speed(0.5))
            .labelled_by(cap.id);
    });
}

/// Draw the Frames Workbench right-side panel. A no-op when the
/// `show_frames_workbench` toggle is off.
pub fn draw_frames_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_frames_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_frames_workbench",
        "Frames",
        |app, ui| {
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
        },
    );
    if close {
        app.show_frames_workbench = false;
    }
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

/// Build the **Frames** result card for the Workbench+Agent bridge — a
/// DATA-ONLY [`crate::WorkspaceProduct`] (`mesh: None`) whose `lines` are the
/// genuine cross-section properties ([`run_frames`]) for the canonical default
/// profile (an IPE 200 I-beam: label, area, perimeter). Registered as the
/// `"frames"` producer in [`crate::products_registry::lookup`]; the tile renders
/// it as a text card, not a 3-D view.
pub(crate) fn frames_product() -> crate::WorkspaceProduct {
    let mut s = FramesWorkbenchState::default();
    run_frames(&mut s);
    crate::WorkspaceProduct {
        title: "Frames".into(),
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
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = FramesWorkbenchState::default();
        s.agent_set("height h", &AgentValue::Float(250.0)).unwrap();
        assert_eq!(s.h, 250.0);
        // `thickness t` and `wall t` both address the `t` field.
        s.agent_set("wall t", &AgentValue::Float(8.0)).unwrap();
        assert_eq!(s.t, 8.0);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into a numeric field) -> Err.
        assert!(s
            .agent_set("height h", &AgentValue::Str("tall".into()))
            .is_err());
        // Non-positive -> Err, field untouched.
        assert!(s.agent_set("height h", &AgentValue::Float(0.0)).is_err());
        assert_eq!(s.h, 250.0, "rejected set leaves field untouched");
    }

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

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::{Node, NodeId, Role};

        // Render with accesskit enabled and read the emitted a11y tree — the
        // same tree a screen reader / AI UI-Automation driver consumes.
        let mut app = ValenxApp::default();
        app.show_frames_workbench = true;
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_frames_workbench(&mut app, ctx);
        });
        let nodes: Vec<(NodeId, Node)> = out
            .platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes;

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // The geometry / loading dimension rows each draw a numeric spin button.
        assert!(
            spin_buttons.len() >= 4,
            "expected the frames numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        // Every DragValue must be associated with a caption (AI-drivable name).
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every frames DragValue must be labelled_by its caption"
        );
    }
}
