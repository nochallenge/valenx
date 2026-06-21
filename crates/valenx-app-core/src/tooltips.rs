//! Tooltip helpers — small extensions to `egui::Response` and a
//! catalogue of standardised tooltip strings the workbench panels
//! share.
//!
//! The polish pass needed hundreds of tooltips spread across the
//! genetics / aero / mesh panels. Rather than scatter free-text
//! strings through every panel module (which makes the wording
//! drift over time), this module exposes:
//!
//! 1. **[`ResponseExt`]** — a trait adding `.tt()` (short for
//!    `.on_hover_text(...)`) and `.tt_with_unit()` so the tooltip
//!    call sites stay terse.
//! 2. **[`run_button_tooltip`]** — a one-line formatter producing
//!    the standard "Runs `{crate}::{fn}`. Shortcut: Ctrl+R / F5"
//!    string. The Genetics workbench's `common::run_button` helper
//!    bakes a similar tooltip into every Run button so the Ctrl+R
//!    hint stays visible everywhere.
//! 3. **[`unit_str`]** — a name → unit lookup the workbenches use
//!    when they want to annotate numeric inputs with their SI unit
//!    in the tooltip ("Speed (m/s)", "Temperature (°C)", …).

use eframe::egui;

/// Convenience trait for adding tooltips to egui responses with less
/// ceremony than the `.on_hover_text(...)` chain. The shorter call
/// site means individual panel functions stay readable when every
/// other control gets a hover annotation.
pub trait ResponseExt {
    /// Short alias for `on_hover_text(text)`.
    fn tt(self, text: impl Into<String>) -> Self;
    /// Tooltip text plus a "Unit: {unit}" line — used by numeric
    /// inputs (DragValue / Slider) so the user can hover to see what
    /// the number means.
    fn tt_with_unit(self, text: impl Into<String>, unit: &str) -> Self;
}

impl ResponseExt for egui::Response {
    fn tt(self, text: impl Into<String>) -> Self {
        self.on_hover_text(text.into())
    }

    fn tt_with_unit(self, text: impl Into<String>, unit: &str) -> Self {
        let text = text.into();
        let combined = format!("{text}\nUnit: {unit}");
        self.on_hover_text(combined)
    }
}

/// Return a tooltip-friendly unit label for a small set of common
/// numeric-input names. Used by the workbench panels' "annotate this
/// DragValue with its unit" calls. Unknown names fall back to an
/// empty string, which the caller treats as "no unit annotation".
pub fn unit_str(name: &str) -> &'static str {
    match name {
        // Aero / wind tunnel.
        "speed" => "m/s",
        "aoa" => "deg",
        "yaw" => "deg",
        "altitude" => "m",
        "temperature" => "°C",
        "pressure" => "Pa",
        "density" => "kg/m³",
        "viscosity" => "Pa·s",
        "reynolds" => "(dimensionless)",
        "mach" => "(dimensionless)",
        "cd" | "cl" | "cm" => "(dimensionless)",
        "force" => "N",
        "moment" => "N·m",
        "area" => "m²",
        "length" | "x" | "y" | "z" => "m",
        // Genetics — temperatures, lengths.
        "tm" => "°C",
        "length_nt" => "nt",
        "length_aa" => "aa",
        "ph" => "(dimensionless)",
        // MD / chemistry.
        "timestep" => "ps",
        "duration_ns" => "ns",
        "kappa" => "1/Å",
        // Time.
        "seconds" => "s",
        "minutes" => "min",
        "iterations" => "iters",
        // Generic.
        "percent" => "%",
        "ratio" => "(0–1)",
        _ => "",
    }
}

/// Format a one-line tooltip for a primary Run button.
///
/// Pattern: `"Run — calls {crate}::{fn}. (Ctrl+R)"`. The Ctrl+R hint
/// shows even on platforms where the user hasn't discovered the
/// shortcut layer yet — the goal is to make the shortcut
/// self-documenting.
pub fn run_button_tooltip(crate_name: &str, fn_name: &str) -> String {
    format!("Runs `{crate_name}::{fn_name}`. Shortcut: Ctrl+R / F5")
}

/// One-line "what this panel does" tooltip used by tab headers.
///
/// The catalogue lives in [`crate::panel_help::short_summary`] —
/// this is a thin re-export so panel-tab call sites have a stable
/// hook even if the summaries get reshuffled.
pub fn panel_tab_tooltip(panel_name: &str) -> String {
    crate::panel_help::short_summary(panel_name).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_str_covers_common_names() {
        assert_eq!(unit_str("speed"), "m/s");
        assert_eq!(unit_str("aoa"), "deg");
        assert_eq!(unit_str("tm"), "°C");
        assert_eq!(unit_str("nonexistent"), "");
    }

    #[test]
    fn run_button_tooltip_includes_shortcut_hint() {
        let s = run_button_tooltip("valenx-bioseq", "translate");
        assert!(s.contains("Ctrl+R"));
        assert!(s.contains("valenx-bioseq"));
        assert!(s.contains("translate"));
    }

    #[test]
    fn panel_tab_tooltip_returns_some_summary() {
        // The panel_help catalogue always returns *some* text — even
        // for an unrecognised panel name the fallback says "no help
        // text available", never empty.
        let s = panel_tab_tooltip("Sequence");
        assert!(!s.is_empty());
        let s = panel_tab_tooltip("UnknownPanel");
        assert!(!s.is_empty());
    }
}
