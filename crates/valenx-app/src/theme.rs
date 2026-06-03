//! Theme system — token-driven egui visual palette.
//!
//! Three variants:
//!
//! 1. **Dark** (default) — the original Valenx CAD-night palette,
//!    sourced from `surface::*` / `text::*` / `accent::*` in
//!    `valenx-design-tokens`.
//! 2. **Light** — a high-luminance counterpart for daytime / projector
//!    use, sourced from `light_surface::*` / `light_text::*` and
//!    accent re-used unchanged (the accents are verified-contrast on
//!    both light and dark surfaces in the contrast-audit test).
//! 3. **High-Contrast** — a WCAG-AAA pure-black + yellow / cyan / green
//!    accent palette intended for low-vision users and bright-sunlight
//!    viewing. Pure black surface, pure white text, saturated yellow
//!    primary accent.
//!
//! The [`ThemeVariant`] enum is the *user-visible* choice; the
//! resolved [`ResolvedTheme`] is what feeds egui (a `Visuals`, a
//! `Style` text size scale, and a few per-call colour accessors
//! the workbenches consume for status badges).
//!
//! ## Why this lives in `valenx-app` and not `valenx-design-tokens`
//!
//! `valenx-design-tokens` is the pure data layer — it doesn't depend
//! on egui. The egui adapter (the `Visuals` builder, the [`apply`]
//! call) lives here because it pulls `eframe::egui`, which we want
//! to keep out of the tokens crate's dependency cone.

use eframe::egui;
use serde::{Deserialize, Serialize};

use valenx_design_tokens::color;

/// Which palette variant the user picked.
///
/// Wire format is kebab-case so settings.json hand-edits read
/// naturally (`"theme-variant": "high-contrast"`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ThemeVariant {
    /// Token-driven dark palette — the default Valenx look.
    #[default]
    Dark,
    /// Token-driven light palette — high-luminance daytime surfaces.
    Light,
    /// WCAG-AAA pure-black + saturated-accent palette for low-vision
    /// users and outdoor / bright-sunlight viewing.
    HighContrast,
}

impl ThemeVariant {
    /// Every variant in display order (matches the radio row in
    /// Settings).
    pub const ALL: [ThemeVariant; 3] = [
        ThemeVariant::Dark,
        ThemeVariant::Light,
        ThemeVariant::HighContrast,
    ];

    /// Short label for the Settings radio.
    pub fn label(self) -> &'static str {
        match self {
            ThemeVariant::Dark => "Dark",
            ThemeVariant::Light => "Light",
            ThemeVariant::HighContrast => "High Contrast",
        }
    }

    /// One-line description for the tooltip.
    pub fn description(self) -> &'static str {
        match self {
            ThemeVariant::Dark => {
                "CAD-night dark palette (default). Verified WCAG-AA contrast on all text tiers."
            }
            ThemeVariant::Light => "High-luminance light palette for daytime / projector use.",
            ThemeVariant::HighContrast => {
                "WCAG-AAA pure-black palette with saturated accents. Designed for low-vision \
                 users and bright-sunlight viewing."
            }
        }
    }
}

/// Resolve a [`ThemeVariant`] into the colour set the egui adapter
/// needs.
#[derive(Clone, Copy, Debug)]
pub struct ResolvedTheme {
    /// Background of the central viewport / window chrome.
    pub surface_0: egui::Color32,
    /// Background of side panels, headers, status bar.
    pub surface_1: egui::Color32,
    /// Background of input fields, dropdowns, button bodies.
    pub surface_2: egui::Color32,
    /// Background of hover / selected rows.
    pub surface_3: egui::Color32,
    /// Body text.
    pub text_1: egui::Color32,
    /// Muted / secondary text.
    pub text_2: egui::Color32,
    /// Disabled / hint text.
    pub text_3: egui::Color32,
    /// Primary accent (selected tab, focus ring, primary button).
    pub accent_primary: egui::Color32,
    /// Success ✓ / green status.
    pub accent_success: egui::Color32,
    /// Warning ⚠ / amber status.
    pub accent_warning: egui::Color32,
    /// Error ✖ / red status.
    pub accent_error: egui::Color32,
    /// Info / informational status.
    pub accent_info: egui::Color32,
}

impl ResolvedTheme {
    /// Build the dark-theme palette from the canonical token set.
    pub fn dark() -> Self {
        Self {
            surface_0: rgb(color::surface::S0),
            surface_1: rgb(color::surface::S1),
            surface_2: rgb(color::surface::S2),
            surface_3: rgb(color::surface::S3),
            text_1: rgb(color::text::T1),
            text_2: rgb(color::text::T2),
            text_3: rgb(color::text::T3),
            accent_primary: rgb(color::accent::PRIMARY),
            accent_success: rgb(color::accent::SUCCESS),
            accent_warning: rgb(color::accent::WARNING),
            accent_error: rgb(color::accent::ERROR),
            accent_info: rgb(color::accent::INFO),
        }
    }

    /// Build the light-theme palette. Surfaces flip to the
    /// `light_surface::*` family, text to `light_text::*`. Accents are
    /// re-used unchanged — verified-contrast on both light and dark
    /// surfaces in `crates/valenx-design-tokens/tests/contrast_audit.rs`.
    pub fn light() -> Self {
        Self {
            surface_0: rgb(color::light_surface::S0),
            surface_1: rgb(color::light_surface::S1),
            surface_2: rgb(color::light_surface::S2),
            surface_3: rgb(color::light_surface::S3),
            text_1: rgb(color::light_text::T1),
            text_2: rgb(color::light_text::T2),
            text_3: rgb(color::light_text::T3),
            accent_primary: rgb(color::accent::PRIMARY),
            accent_success: rgb(color::accent::SUCCESS),
            accent_warning: rgb(color::accent::WARNING),
            accent_error: rgb(color::accent::ERROR),
            accent_info: rgb(color::accent::INFO),
        }
    }

    /// Build the high-contrast palette: pure black surfaces, pure
    /// white text, saturated yellow / cyan accents. Each accent
    /// achieves > 7:1 contrast on `hc_surface::S1` (pure black),
    /// well above WCAG-AAA for normal text.
    pub fn high_contrast() -> Self {
        Self {
            surface_0: rgb(color::hc_surface::S0),
            surface_1: rgb(color::hc_surface::S1),
            surface_2: rgb(color::hc_surface::S2),
            surface_3: rgb(color::hc_surface::S3),
            text_1: rgb(color::hc_text::T1),
            text_2: rgb(color::hc_text::T2),
            text_3: rgb(color::hc_text::T3),
            accent_primary: rgb(color::hc_accent::PRIMARY),
            accent_success: rgb(color::hc_accent::SUCCESS),
            accent_warning: rgb(color::hc_accent::WARNING),
            accent_error: rgb(color::hc_accent::ERROR),
            accent_info: rgb(color::hc_accent::INFO),
        }
    }

    /// Resolve a variant.
    pub fn from_variant(v: ThemeVariant) -> Self {
        match v {
            ThemeVariant::Dark => Self::dark(),
            ThemeVariant::Light => Self::light(),
            ThemeVariant::HighContrast => Self::high_contrast(),
        }
    }
}

/// Apply the resolved theme + font scale to an egui context.
///
/// Compose: builds a `Visuals` from the resolved palette and a
/// `Style` whose text-style sizes are scaled by `font_scale` (1.0 =
/// default; range 0.75–2.0 per the Settings dialog clamp).
pub fn apply(ctx: &egui::Context, variant: ThemeVariant, font_scale: f32) {
    let resolved = ResolvedTheme::from_variant(variant);

    // Start from egui's dark / light defaults so we inherit sensible
    // values for the dozens of secondary colours we don't override
    // (selection underline, scroll handle, etc.).
    let mut visuals = if matches!(variant, ThemeVariant::Light) {
        egui::Visuals::light()
    } else {
        egui::Visuals::dark()
    };

    // Window / panel / central backdrop.
    visuals.window_fill = resolved.surface_1;
    visuals.panel_fill = resolved.surface_1;
    visuals.extreme_bg_color = resolved.surface_0;
    visuals.faint_bg_color = resolved.surface_2;

    // Widget body fills — inactive / active / hovered states use the
    // surface tiers so the entire UI keeps a coherent depth scale.
    visuals.widgets.noninteractive.bg_fill = resolved.surface_1;
    visuals.widgets.inactive.bg_fill = resolved.surface_2;
    visuals.widgets.hovered.bg_fill = resolved.surface_3;
    visuals.widgets.active.bg_fill = resolved.surface_3;
    visuals.widgets.open.bg_fill = resolved.surface_3;

    // Body text + secondary text. egui pulls override_text_color when
    // building per-widget text colours, so a single hook covers most
    // labels / buttons / dropdowns.
    visuals.override_text_color = Some(resolved.text_1);
    visuals.widgets.noninteractive.fg_stroke.color = resolved.text_1;
    visuals.widgets.inactive.fg_stroke.color = resolved.text_1;
    visuals.widgets.hovered.fg_stroke.color = resolved.text_1;
    visuals.widgets.active.fg_stroke.color = resolved.text_1;
    visuals.widgets.open.fg_stroke.color = resolved.text_1;

    // Accent — egui's hyperlink_color is the canonical "this is the
    // primary accent" hook (used by selectable_value highlights, focus
    // rings, and text-link colour).
    visuals.hyperlink_color = resolved.accent_primary;
    visuals.selection.bg_fill = resolved.accent_primary.linear_multiply(0.35);
    visuals.selection.stroke.color = resolved.accent_primary;

    // High-contrast variant: thicken stroke widths for outlines and
    // bump the focus stroke so keyboard-focused controls are visible
    // from across the room.
    if matches!(variant, ThemeVariant::HighContrast) {
        visuals.widgets.noninteractive.bg_stroke.width = 1.5;
        visuals.widgets.inactive.bg_stroke.width = 1.5;
        visuals.widgets.hovered.bg_stroke.width = 2.0;
        visuals.widgets.active.bg_stroke.width = 2.0;
        visuals.widgets.noninteractive.bg_stroke.color = resolved.text_1;
        visuals.widgets.inactive.bg_stroke.color = resolved.text_1;
        visuals.widgets.hovered.bg_stroke.color = resolved.accent_primary;
        visuals.widgets.active.bg_stroke.color = resolved.accent_primary;
        // Pure-white text on pure-black is the AAA gold standard —
        // we don't want egui to muddy it with a tint.
        visuals.override_text_color = Some(resolved.text_1);
    }

    ctx.set_visuals(visuals);

    // Font-scale slider — clamps in the Settings dialog at 0.75–2.0
    // so values outside that band can't reach this branch via
    // normal UI. The defensive clamp here covers a hand-edited
    // settings.json file.
    let clamped = font_scale.clamp(0.5, 3.0);
    let mut style: egui::Style = (*ctx.style()).clone();
    for (_, font_id) in style.text_styles.iter_mut() {
        // The default text-style sizes come from the egui default
        // style — scale them in place. clamped is multiplicative so
        // 1.0 leaves the defaults untouched.
        font_id.size *= clamped;
    }

    // --- Layout polish: a calmer, more deliberate spacing + subtle
    // rounding than egui's tight defaults. Applied centrally here so every
    // panel and control inherits it — the single biggest "commercial feel"
    // lever without touching individual panels. Tuned to add breathing
    // room without bloating the dense toolbox into excessive scrolling.
    style.spacing.item_spacing = egui::vec2(8.0, 5.0);
    style.spacing.button_padding = egui::vec2(8.0, 4.0);
    style.spacing.menu_margin = egui::Margin::same(6.0);
    style.spacing.indent = 18.0;
    style.spacing.interact_size.y = 22.0;
    let widget_rounding = egui::Rounding::same(4.0);
    style.visuals.widgets.noninteractive.rounding = widget_rounding;
    style.visuals.widgets.inactive.rounding = widget_rounding;
    style.visuals.widgets.hovered.rounding = widget_rounding;
    style.visuals.widgets.active.rounding = widget_rounding;
    style.visuals.widgets.open.rounding = widget_rounding;
    style.visuals.window_rounding = egui::Rounding::same(6.0);
    style.visuals.menu_rounding = egui::Rounding::same(6.0);

    ctx.set_style(style);
}

/// Convert a `0xRRGGBB` token const to an egui colour.
fn rgb(hex: u32) -> egui::Color32 {
    egui::Color32::from_rgb(
        ((hex >> 16) & 0xFF) as u8,
        ((hex >> 8) & 0xFF) as u8,
        (hex & 0xFF) as u8,
    )
}

/// Clamp a font scale to the supported UI range.
///
/// The range matches the Settings slider: 0.75–2.0 covers everything
/// from "fits more on screen" to "AAA-large body text". Values
/// outside the range either pack the UI too tight to be clickable
/// (< 0.75) or push popovers off-screen (> 2.0).
pub fn clamp_font_scale(s: f32) -> f32 {
    s.clamp(0.75, 2.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variant_all_lists_every_case() {
        assert_eq!(ThemeVariant::ALL.len(), 3);
        for v in ThemeVariant::ALL {
            assert!(!v.label().is_empty());
            assert!(!v.description().is_empty());
        }
    }

    #[test]
    fn default_variant_is_dark() {
        assert_eq!(ThemeVariant::default(), ThemeVariant::Dark);
    }

    #[test]
    fn variant_serializes_kebab_case() {
        let v = ThemeVariant::HighContrast;
        let s = serde_json::to_string(&v).unwrap();
        assert_eq!(s, "\"high-contrast\"");
        let back: ThemeVariant = serde_json::from_str(&s).unwrap();
        assert_eq!(back, v);
    }

    #[test]
    fn resolve_each_variant_produces_distinct_palette() {
        let d = ResolvedTheme::dark();
        let l = ResolvedTheme::light();
        let hc = ResolvedTheme::high_contrast();
        // The three palettes must not all collapse to the same colour
        // — that would indicate the token wiring is broken.
        assert_ne!(d.surface_0, l.surface_0);
        assert_ne!(d.surface_0, hc.surface_0);
        assert_ne!(l.surface_0, hc.surface_0);
        // Pure-black is the high-contrast surface_0 marker.
        assert_eq!(hc.surface_0, egui::Color32::BLACK);
    }

    #[test]
    fn high_contrast_text_is_pure_white() {
        let hc = ResolvedTheme::high_contrast();
        assert_eq!(hc.text_1, egui::Color32::WHITE);
    }

    #[test]
    fn clamp_font_scale_obeys_range() {
        assert_eq!(clamp_font_scale(0.5), 0.75);
        assert_eq!(clamp_font_scale(1.0), 1.0);
        assert_eq!(clamp_font_scale(3.0), 2.0);
        assert_eq!(clamp_font_scale(1.5), 1.5);
    }

    #[test]
    fn apply_does_not_panic_on_any_variant() {
        let ctx = egui::Context::default();
        for v in ThemeVariant::ALL {
            apply(&ctx, v, 1.0);
        }
    }

    #[test]
    fn apply_respects_font_scale() {
        // After applying scale 2.0, every text-style size should be
        // about twice the baseline.
        let ctx = egui::Context::default();
        let baseline_sizes: Vec<f32> =
            ctx.style().text_styles.values().map(|f| f.size).collect();
        apply(&ctx, ThemeVariant::Dark, 2.0);
        let scaled_sizes: Vec<f32> =
            ctx.style().text_styles.values().map(|f| f.size).collect();
        assert_eq!(baseline_sizes.len(), scaled_sizes.len());
        for (b, s) in baseline_sizes.iter().zip(scaled_sizes.iter()) {
            assert!((s / b - 2.0).abs() < 0.01, "expected 2x scale: {b} → {s}");
        }
    }
}
