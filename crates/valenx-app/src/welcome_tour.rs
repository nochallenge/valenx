//! First-run welcome tour — a 3-step orientation popup.
//!
//! Independent from [`crate::first_run`] (which probes the adapter
//! environment): this is a short visual walkthrough of the three
//! workbenches that runs on top of the wizard so a brand-new user
//! lands on a screen that explains what Valenx actually is.
//!
//! The tour is suppressed once dismissed via the
//! `Settings.welcome_tour_completed` flag (settings.json). Re-openable
//! from the Help menu and from the command palette.

use eframe::egui;

/// Tour state — which step the user is on (0..=2) plus a "finished"
/// flag the host inspects after [`render`] to decide whether to
/// persist the welcome-tour-completed bit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct TourState {
    /// 0-indexed step. The tour runs through steps 0, 1, 2; on a
    /// Next click from step 2 it sets `finished = true`.
    pub step: u8,
    /// `true` once the user clicked "Get started" on the final step
    /// or "Skip" on any step.
    pub finished: bool,
}

impl TourState {
    /// Begin a fresh tour at step 0.
    pub fn new() -> Self {
        Self::default()
    }
}

/// Number of tour steps.
pub const STEP_COUNT: u8 = 3;

/// One-line title for the step.
pub fn step_title(step: u8) -> &'static str {
    match step {
        0 => "Welcome to Valenx",
        1 => "Four workbenches",
        2 => "Shortcuts to know",
        _ => "Ready",
    }
}

/// Multiline body for the step.
pub fn step_body(step: u8) -> &'static str {
    match step {
        0 => {
            "\
Valenx is a native open-source desktop suite that unifies CAD, CFD,
FEA, electromagnetics, chemistry, molecular dynamics, computational
biology, and external-aerodynamics into a single Rust shell.

No browser, no subscription, no vendor lock-in.

This 3-step tour will orient you in under a minute.
"
        }
        1 => {
            "\
The right-hand side of the window has four workbenches:

  Ctrl+1 — Mesh Toolbox (CAD: Part, Draft, TechDraw, Assembly,
           Surface, CAM, Architecture, Sketcher, Part Design …).

  Ctrl+2 — Genetics workbench (15 panels covering sequence editing,
           alignment, phylogenetics, RNA structure, MD,
           cheminformatics, docking, gene editing …).

  Ctrl+3 — Wind Tunnel workbench (3-D external-aerodynamics CFD —
           drag / lift / moment, Cp + velocity overlays, AoA polar).

  Ctrl+4 — Astro / Launch workbench (launch ascent to orbit — Δv
           budget, max-Q, staging + flight profile — and Hohmann /
           hoverslam / rendezvous / launch-azimuth planners).

Each workbench is independent; egui docks them side by side when
several are open.
"
        }
        2 => {
            "\
A handful of shortcuts unlock the rest of the UI:

  Ctrl+P — Command palette (fuzzy-search every action).
  Ctrl+R — Run the active panel's primary action.
  Ctrl+Z / Ctrl+Y — Undo / redo edits.
  F1 — Contextual help for the focused panel.
  ?  — Toggle this cheat-sheet overlay any time.

You can change the colour palette (Dark / Light / High-Contrast)
and the font scale from Settings → Appearance.

Have fun.
"
        }
        _ => "",
    }
}

/// Render the welcome-tour popup. Mutates `state.step` / `.finished`
/// in place. Returns `true` if the tour was dismissed this frame
/// (the host saves the completed bit + drops the open flag).
pub fn render(ctx: &egui::Context, open: &mut bool, state: &mut TourState) -> bool {
    let mut dismissed_this_frame = false;
    egui::Window::new(format!(
        "{} ({}/{})",
        step_title(state.step),
        state.step + 1,
        STEP_COUNT
    ))
    .open(open)
    .collapsible(false)
    .resizable(false)
    .default_width(440.0)
    .show(ctx, |ui| {
        ui.add_space(4.0);
        ui.label(step_body(state.step));
        ui.add_space(8.0);
        ui.separator();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(state.step > 0, egui::Button::new("Back"))
                .on_hover_text("Previous step")
                .clicked()
            {
                state.step = state.step.saturating_sub(1);
            }

            let on_last = state.step + 1 == STEP_COUNT;
            let primary_label = if on_last { "Get started" } else { "Next" };
            if ui
                .button(primary_label)
                .on_hover_text(if on_last {
                    "Finish the tour and remember the choice."
                } else {
                    "Continue to the next step."
                })
                .clicked()
            {
                if on_last {
                    state.finished = true;
                    dismissed_this_frame = true;
                } else {
                    state.step += 1;
                }
            }

            if ui
                .button("Skip")
                .on_hover_text("Dismiss the tour. Re-openable from Help → Welcome tour.")
                .clicked()
            {
                state.finished = true;
                dismissed_this_frame = true;
            }
        });
    });

    if dismissed_this_frame {
        *open = false;
    }
    dismissed_this_frame
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fresh_state_starts_at_step_zero() {
        let s = TourState::new();
        assert_eq!(s.step, 0);
        assert!(!s.finished);
    }

    #[test]
    fn step_metadata_exists_for_every_step() {
        for step in 0..STEP_COUNT {
            assert!(!step_title(step).is_empty());
            assert!(!step_body(step).is_empty());
        }
    }

    #[test]
    fn step_count_matches_metadata() {
        // The tour ships with STEP_COUNT entries — bumping the count
        // without filling in title / body would render a blank popup.
        for step in 0..STEP_COUNT {
            let body = step_body(step);
            assert!(body.lines().count() >= 2);
        }
    }

    #[test]
    fn render_does_not_panic_at_any_step() {
        let ctx = egui::Context::default();
        let mut open = true;
        let mut state = TourState::new();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            for step in 0..STEP_COUNT {
                state.step = step;
                let _ = render(ctx, &mut open, &mut state);
            }
        });
    }
}
