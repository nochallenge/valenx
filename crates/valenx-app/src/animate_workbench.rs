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

/// Draw the animation workbench (a no-op unless toggled on via
/// View → Animation).
pub fn draw_animate_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_animate_workbench {
        return;
    }
    egui::SidePanel::right("valenx_animate_workbench")
        .resizable(true)
        .default_width(380.0)
        .width_range(320.0..=720.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(ui, "Animation", "keyframe timeline · valenx-animate") {
                app.show_animate_workbench = false;
            }
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
        });
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
