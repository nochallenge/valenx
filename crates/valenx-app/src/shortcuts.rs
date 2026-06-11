//! Keyboard-shortcut layer.
//!
//! The legacy `update.rs` had inline `i.modifiers.ctrl && i.key_pressed(...)`
//! checks scattered through the input handler. The polish pass needed
//! ~ a dozen new shortcuts (workbench switching, Run, Save, Undo,
//! Redo, F1, ?, etc.) — adding them inline would have made the
//! handler unreadable. This module hoists the entire shortcut
//! catalogue into one place so:
//!
//! 1. The canonical list of bindings is **a single `enum`** the
//!    keyboard-cheat-sheet overlay can iterate to render itself —
//!    one source of truth, no drift between "what is bound" and
//!    "what the help overlay claims is bound".
//! 2. A test can verify that every variant has a distinct binding
//!    (no accidental collisions) and a non-empty description.
//! 3. Callers consume a typed [`ShortcutAction`] rather than
//!    matching on raw `egui::Key` values inside the giant `update`
//!    fn.
//!
//! ## What this module does NOT do
//!
//! It doesn't *dispatch* the actions — that's the host's job. The
//! input-handler in `update.rs` calls [`poll_shortcut`] and gets back
//! `Vec<ShortcutAction>`; then a `match` in the caller routes each
//! action to the right `ValenxApp` method (which may need `&mut self`
//! that the input-closure can't hold).

use eframe::egui;

use crate::genetics_workbench::GeneticsPanel;

/// A high-level action a keyboard shortcut maps to. The host matches
/// on these inside `update` to call the appropriate `ValenxApp`
/// method.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ShortcutAction {
    // ---- Workbench switching ----------------------------------------
    /// Show the right-side Mesh Toolbox panel.
    ShowMeshToolbox,
    /// Show the right-side Genetics Workbench panel.
    ShowGeneticsWorkbench,
    /// Show the right-side Aero / Wind Tunnel workbench panel.
    ShowAeroWorkbench,
    /// Show the right-side Astro / Launch workbench panel.
    ShowAstroWorkbench,

    // ---- File / project ---------------------------------------------
    /// Open the "File → New Project…" modal so the user can scaffold
    /// a fresh project from a template without dropping to a CLI.
    NewProject,
    /// Open the native folder picker for an existing `.valenx`
    /// directory. Mirrors File → Open project.
    OpenProject,

    // ---- Common actions ---------------------------------------------
    /// Trigger the Run / Compute action for the currently visible
    /// workbench's primary action. The host inspects which workbench
    /// is up and calls the right method.
    RunPrimary,
    /// Save the current state.
    SavePrimary,
    /// Undo the last reversible edit in the active panel.
    Undo,
    /// Redo the previously undone edit.
    Redo,
    /// Cancel the in-flight long-running operation (run / sweep /
    /// aero solve).
    Cancel,

    // ---- Help / palette ---------------------------------------------
    /// Open the command palette overlay.
    OpenCommandPalette,
    /// Toggle the keyboard-shortcut cheat-sheet overlay.
    ToggleKeyboardHelp,
    /// Open the contextual help popup for the currently focused
    /// panel.
    OpenContextHelp,
}

impl ShortcutAction {
    /// Every action in display order. Used by the keyboard-cheat-sheet
    /// overlay to render itself.
    pub const ALL: [ShortcutAction; 14] = [
        ShortcutAction::OpenCommandPalette,
        ShortcutAction::OpenContextHelp,
        ShortcutAction::ToggleKeyboardHelp,
        ShortcutAction::NewProject,
        ShortcutAction::OpenProject,
        ShortcutAction::ShowMeshToolbox,
        ShortcutAction::ShowGeneticsWorkbench,
        ShortcutAction::ShowAeroWorkbench,
        ShortcutAction::ShowAstroWorkbench,
        ShortcutAction::RunPrimary,
        ShortcutAction::SavePrimary,
        ShortcutAction::Undo,
        ShortcutAction::Redo,
        ShortcutAction::Cancel,
    ];

    /// Short label for the cheat-sheet (e.g. "Open command palette").
    pub fn label(self) -> &'static str {
        match self {
            ShortcutAction::OpenCommandPalette => "Open command palette",
            ShortcutAction::OpenContextHelp => "Show this panel's help",
            ShortcutAction::ToggleKeyboardHelp => "Toggle keyboard cheat-sheet",
            ShortcutAction::NewProject => "New project",
            ShortcutAction::OpenProject => "Open project",
            ShortcutAction::ShowMeshToolbox => "Show Mesh Toolbox (workbench 1)",
            ShortcutAction::ShowGeneticsWorkbench => "Show Genetics workbench (workbench 2)",
            ShortcutAction::ShowAeroWorkbench => "Show Wind Tunnel workbench (workbench 3)",
            ShortcutAction::ShowAstroWorkbench => "Show Astro / Launch workbench (workbench 4)",
            ShortcutAction::RunPrimary => "Run / Compute the active panel's primary action",
            ShortcutAction::SavePrimary => "Save (project / settings / panel state)",
            ShortcutAction::Undo => "Undo the last reversible edit",
            ShortcutAction::Redo => "Redo the undone edit",
            ShortcutAction::Cancel => "Cancel the in-flight run / sweep",
        }
    }

    /// Display string for the binding, e.g. `"Ctrl+1"`, `"F1"`,
    /// `"?"`. Used by the cheat-sheet.
    pub fn binding(self) -> &'static str {
        match self {
            ShortcutAction::OpenCommandPalette => "Ctrl+P",
            ShortcutAction::OpenContextHelp => "F1",
            ShortcutAction::ToggleKeyboardHelp => "?",
            ShortcutAction::NewProject => "Ctrl+N",
            ShortcutAction::OpenProject => "Ctrl+O",
            ShortcutAction::ShowMeshToolbox => "Ctrl+1",
            ShortcutAction::ShowGeneticsWorkbench => "Ctrl+2",
            ShortcutAction::ShowAeroWorkbench => "Ctrl+3",
            ShortcutAction::ShowAstroWorkbench => "Ctrl+4",
            ShortcutAction::RunPrimary => "Ctrl+R",
            ShortcutAction::SavePrimary => "Ctrl+S",
            ShortcutAction::Undo => "Ctrl+Z",
            ShortcutAction::Redo => "Ctrl+Y",
            ShortcutAction::Cancel => "Esc",
        }
    }

    /// Per-action one-line description for the help-overlay's body
    /// column (alongside the label + binding).
    pub fn description(self) -> &'static str {
        match self {
            ShortcutAction::OpenCommandPalette => {
                "Fuzzy-search every action across all workbenches."
            }
            ShortcutAction::OpenContextHelp => {
                "Pop up help text for the currently focused panel."
            }
            ShortcutAction::ToggleKeyboardHelp => {
                "Open / close the keyboard-shortcut cheat-sheet."
            }
            ShortcutAction::NewProject => {
                "Open the New Project modal — pick a template and a location, scaffold the .valenx directory."
            }
            ShortcutAction::OpenProject => {
                "Open the native folder picker for an existing .valenx project."
            }
            ShortcutAction::ShowMeshToolbox => {
                "Switch to the Mesh Toolbox (Part / Draft / Assembly / CAM / etc.)."
            }
            ShortcutAction::ShowGeneticsWorkbench => {
                "Switch to the Genetics workbench (15 computational-biology panels)."
            }
            ShortcutAction::ShowAeroWorkbench => {
                "Switch to the Wind Tunnel workbench (3-D external-aero CFD)."
            }
            ShortcutAction::ShowAstroWorkbench => {
                "Switch to the Astro / Launch workbench (launch ascent sim + mission planners)."
            }
            ShortcutAction::RunPrimary => {
                "Triggers the active panel's main button (Translate / Fold / Solve …)."
            }
            ShortcutAction::SavePrimary => "Save the project / settings (where applicable).",
            ShortcutAction::Undo => "Reverse the most recent edit in the active panel.",
            ShortcutAction::Redo => "Reapply the undone edit.",
            ShortcutAction::Cancel => {
                "Stop a long-running solve / sweep / run. Closes overlays first."
            }
        }
    }
}

/// Resolve the keyboard input for this frame into the set of
/// shortcut-actions it fired. Called from inside an
/// `ctx.input(|i| { ... })` closure so the host doesn't have to
/// open / close the input lock twice.
///
/// Avoids returning the action twice if multiple bindings overlap
/// — every binding's key check returns once.
pub fn poll_shortcut(i: &egui::InputState, palette_open: bool) -> Vec<ShortcutAction> {
    let mut out = Vec::new();
    let ctrl = i.modifiers.ctrl || i.modifiers.command;

    // Palette has priority; if it's open we don't want F1 / Ctrl+R
    // to leak through (the palette swallows arrows + enter + esc on
    // its own). Esc + the keyboard-help "?" toggle still get through
    // so the user can dismiss overlays while the palette is up.
    if !palette_open {
        if ctrl && i.key_pressed(egui::Key::P) {
            out.push(ShortcutAction::OpenCommandPalette);
        }
        if ctrl && i.key_pressed(egui::Key::N) {
            out.push(ShortcutAction::NewProject);
        }
        if ctrl && i.key_pressed(egui::Key::O) {
            out.push(ShortcutAction::OpenProject);
        }
        if ctrl && i.key_pressed(egui::Key::Num1) {
            out.push(ShortcutAction::ShowMeshToolbox);
        }
        if ctrl && i.key_pressed(egui::Key::Num2) {
            out.push(ShortcutAction::ShowGeneticsWorkbench);
        }
        if ctrl && i.key_pressed(egui::Key::Num3) {
            out.push(ShortcutAction::ShowAeroWorkbench);
        }
        if ctrl && i.key_pressed(egui::Key::Num4) {
            out.push(ShortcutAction::ShowAstroWorkbench);
        }
        if ctrl && i.key_pressed(egui::Key::R) {
            out.push(ShortcutAction::RunPrimary);
        }
        if ctrl && i.key_pressed(egui::Key::S) {
            out.push(ShortcutAction::SavePrimary);
        }
        if ctrl && i.key_pressed(egui::Key::Z) && !i.modifiers.shift {
            out.push(ShortcutAction::Undo);
        }
        // Two ways to redo — Ctrl+Y (Windows-y) and Ctrl+Shift+Z
        // (Mac / Linux-y). Both fire the same action.
        if ctrl && i.key_pressed(egui::Key::Y) {
            out.push(ShortcutAction::Redo);
        }
        if ctrl && i.modifiers.shift && i.key_pressed(egui::Key::Z) {
            out.push(ShortcutAction::Redo);
        }
        if i.key_pressed(egui::Key::F1) {
            out.push(ShortcutAction::OpenContextHelp);
        }
    }

    // The "?" toggle (Shift+/) reaches the host even when the palette
    // is open so the user can pull up the cheat-sheet from any state.
    // egui delivers it via the Slash key + shift modifier (or as a
    // `Key::Questionmark` on some layouts — we match either).
    if (i.modifiers.shift && i.key_pressed(egui::Key::Slash))
        || i.key_pressed(egui::Key::Questionmark)
    {
        out.push(ShortcutAction::ToggleKeyboardHelp);
    }

    if i.key_pressed(egui::Key::Escape) {
        out.push(ShortcutAction::Cancel);
    }

    out
}

/// Map a Ctrl+digit workbench shortcut to the genetics panel it
/// targets — only used when the genetics workbench is the active
/// one. Returns `None` for any digit outside the panel range.
pub fn genetics_panel_for_digit(digit: u8) -> Option<GeneticsPanel> {
    let panels = GeneticsPanel::ALL;
    let idx = (digit as usize).checked_sub(1)?;
    panels.get(idx).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_actions_have_unique_bindings() {
        // No two actions can share the same binding string. This
        // catches accidental copy-paste collisions when adding new
        // shortcuts.
        let mut seen = std::collections::HashSet::new();
        for a in ShortcutAction::ALL {
            assert!(
                seen.insert(a.binding()),
                "{a:?} shares its binding with another action: {}",
                a.binding()
            );
        }
    }

    #[test]
    fn all_actions_have_non_empty_metadata() {
        for a in ShortcutAction::ALL {
            assert!(!a.label().is_empty(), "{a:?} has an empty label");
            assert!(!a.binding().is_empty(), "{a:?} has an empty binding");
            assert!(!a.description().is_empty(), "{a:?} has an empty description");
        }
    }

    #[test]
    fn genetics_panel_for_digit_clamps_to_range() {
        // Track the panel count via ALL so adding a panel doesn't break
        // this test (it did once, when the count was hard-coded at 14).
        let n = GeneticsPanel::ALL.len() as u8;
        assert_eq!(genetics_panel_for_digit(1), Some(GeneticsPanel::Sequence));
        assert!(
            genetics_panel_for_digit(n).is_some(),
            "the last panel (digit {n}) must map"
        );
        assert!(
            genetics_panel_for_digit(n + 1).is_none(),
            "one past the end (digit {}) must not map",
            n + 1
        );
        assert!(genetics_panel_for_digit(0).is_none());
    }

    #[test]
    fn poll_shortcut_returns_empty_on_idle_input() {
        // An empty InputState (no keys pressed) must produce no
        // shortcut actions — guards against accidental fall-through
        // in the polling routine.
        let i = egui::InputState::default();
        let out = poll_shortcut(&i, false);
        assert!(out.is_empty());
    }
}
