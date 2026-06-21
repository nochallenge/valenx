//! Animation workbench — keyframe timeline on `valenx-animate`.
//!
//! A right-side panel that builds a small joint **keyframe animation** (one
//! revolute joint sweeping 0 → π over two seconds, with a selectable easing
//! curve), plots the sampled joint value across time, and exposes a playhead
//! slider with a live readout. The animation sampling is headless-testable.

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};
use std::f64::consts::PI;

use valenx_animate::{Animation, Keyframe, TweenMode};

use crate::ValenxApp;

/// Persistent state for the animation workbench.
pub struct AnimateWorkbenchState {
    tween: TweenMode,
    anim: Animation,
    playhead: f64,
}

impl Default for AnimateWorkbenchState {
    fn default() -> Self {
        let tween = TweenMode::EaseInOut;
        Self {
            anim: demo_animation(tween),
            tween,
            playhead: 0.0,
        }
    }
}

/// A two-keyframe demo: joint 0 sweeps 0 → π over 2 s with the given easing.
fn demo_animation(tween: TweenMode) -> Animation {
    let mut a = Animation::new();
    let _ = a.push(Keyframe::at(0.0).with_joint(0, 0.0));
    let _ = a.push(Keyframe::at(2.0).with_joint(0, PI).with_tween(tween));
    a
}

/// Sample joint 0's value (rad) at time `t`.
fn joint0(anim: &Animation, t: f64) -> f64 {
    anim.sample(t)
        .ok()
        .and_then(|m| m.get(&0).copied())
        .unwrap_or(0.0)
}

/// Build the agent-bridge **`animate` product** — a DATA-ONLY *text card*
/// summarising the keyframe timeline (`mesh: None`, populated `lines`).
///
/// Rendering a single *posed frame* as an image is impractical here: this
/// workbench animates a bare revolute-joint *angle* over time (the
/// `valenx-animate` `Animation` samples joint values, not a posed renderable
/// mesh), so there is no body to rasterise into a frame image. The honest
/// wiring is therefore the same mesh-less text-card path the DATA-ONLY
/// workbenches (`cfd` / `astro` / `car`) use: the card reports the genuine
/// timeline facts — keyframe count, duration, easing, and the sampled
/// start / mid / end joint values — computed from the canonical demo
/// animation. (Distinct from the `render` product, which *does* produce a real
/// raster image.)
pub(crate) fn animate_product() -> crate::WorkspaceProduct {
    let tween = TweenMode::EaseInOut;
    let anim = demo_animation(tween);
    let dur = anim.duration();
    let n = anim.keyframes.len();
    let start = joint0(&anim, 0.0);
    let mid = joint0(&anim, dur * 0.5);
    let end = joint0(&anim, dur);
    crate::WorkspaceProduct {
        title: "Keyframe animation".into(),
        lines: vec![
            format!("{n} keyframes · {dur:.2} s · {} easing", tween.label()),
            "joint 0 sweep (rad):".into(),
            format!("  t=0.00 s  →  {start:.3}"),
            format!("  t={:.2} s  →  {mid:.3}", dur * 0.5),
            format!("  t={dur:.2} s  →  {end:.3}"),
            "valenx-animate · keyframe timeline".into(),
        ],
        mesh: None,
        vertex_colors: None,
        camera: valenx_viz::OrbitCamera::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

/// Draw the animation workbench (a no-op unless toggled on via
/// View → Animation).
pub fn draw_animate_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_animate_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_animate_workbench",
        "Animation",
        |app, ui| {
            ui.label(
                egui::RichText::new("keyframe timeline · valenx-animate")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.animate;
            let mut tw = s.tween;
            egui::ComboBox::from_label("easing")
                .selected_text(tw.label())
                .show_ui(ui, |ui| {
                    for t in [
                        TweenMode::Linear,
                        TweenMode::EaseIn,
                        TweenMode::EaseOut,
                        TweenMode::EaseInOut,
                        TweenMode::Cubic,
                    ] {
                        ui.selectable_value(&mut tw, t, t.label());
                    }
                });
            if tw != s.tween {
                s.tween = tw;
                s.anim = demo_animation(tw);
            }
            let dur = s.anim.duration().max(0.001);
            ui.add(egui::Slider::new(&mut s.playhead, 0.0..=dur).text("time (s)"));
            let v = joint0(&s.anim, s.playhead);
            ui.label(
                egui::RichText::new(format!("joint 0 = {v:.3} rad @ t = {:.2} s", s.playhead))
                    .monospace()
                    .small(),
            );
            ui.separator();
            ui.label(egui::RichText::new("Joint value over time").strong());
            Plot::new("animate_curve").height(160.0).show(ui, |pui| {
                let n = 120;
                let pts: Vec<[f64; 2]> = (0..=n)
                    .map(|i| {
                        let t = dur * i as f64 / n as f64;
                        [t, joint0(&s.anim, t)]
                    })
                    .collect();
                pui.line(Line::new(PlotPoints::from(pts)).name("joint 0 (rad)"));
            });
        },
    );
    if close {
        app.show_animate_workbench = false;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_animation_samples_its_endpoints() {
        let a = demo_animation(TweenMode::Linear);
        assert!(joint0(&a, 0.0).abs() < 1e-9, "starts at 0");
        assert!((joint0(&a, 2.0) - PI).abs() < 1e-6, "ends at π");
    }

    #[test]
    fn midpoint_is_between_endpoints() {
        let a = demo_animation(TweenMode::Linear);
        let mid = joint0(&a, 1.0);
        assert!(mid > 0.0 && mid < PI, "midpoint between endpoints: {mid}");
    }
}
