//! App preferences. Persists to `<state_dir>/settings.json` so user
//! choices survive app restarts.
//!
//! The set is deliberately tiny — just the knobs that affect how the
//! current UI looks or behaves. Physics-related knobs live on the
//! case, not here.

use eframe::egui;
use serde::{Deserialize, Serialize};

use crate::theme::{self, ThemeVariant};
use crate::viewport::ShadingMode;

/// Which colour scheme the app uses. `Auto` matches the OS where
/// supported, else falls back to `Dark` (Valenx is a CAD app; the
/// default tone is deliberately dark so viewport contrast is high).
///
/// Wire format is lowercase (`"auto"` / `"dark"` / `"light"`) — both
/// to match the JSON convention and to be friendly to hand-edited
/// `settings.json` files.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Theme {
    #[default]
    Auto,
    Dark,
    Light,
}

/// Which scale the residual chart uses. Log is the right answer for
/// CFD residuals (which cover 6+ orders of magnitude in one run);
/// linear is occasionally useful for quick diagnostics.
///
/// Wire format is lowercase (`"log10"` / `"linear"`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ResidualScale {
    #[default]
    Log10,
    Linear,
}

/// User-tunable app preferences.
///
/// `Default` is implemented manually below (instead of derived) so
/// `convergence_target` defaults to `Some(1e-4)` everywhere — the
/// derived `Default` would use `Option::None` and the on-disk
/// serde-default would only fire during deserialisation, not for
/// in-memory `Settings::default()` calls.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Settings {
    #[serde(default)]
    pub theme: Theme,
    #[serde(default)]
    pub default_shading: ShadingMode,
    #[serde(default)]
    pub residual_scale: ResidualScale,
    /// Convergence target overlay on the residual chart. When
    /// `Some(threshold)`, a horizontal reference line is drawn at
    /// log10(threshold) (or `threshold` for linear mode); residual
    /// curves crossing it visually flag converged fields. Common
    /// CFD defaults: `1e-4` for engineering accuracy, `1e-6` for
    /// research-grade.
    #[serde(default = "default_convergence_target")]
    pub convergence_target: Option<f64>,
    /// When true, the adapter probe runs every time the Settings
    /// dialog closes. Useful if the user just installed a tool and
    /// wants the registry to pick it up without restarting.
    #[serde(default)]
    pub reprobe_on_close: bool,
    /// Crash-report upload opt-in. **Defaults to `false`** — the
    /// crash reporter (`valenx-crash-reporter`) always writes
    /// reports to `<state_dir>/crashes/` so users can review them
    /// locally; this flag controls whether the next launch's
    /// "found N unsent crash reports — submit them?" prompt
    /// auto-sends rather than asks.
    ///
    /// Network egress is the user's affirmative choice; we never
    /// upload silently.
    #[serde(default)]
    pub crash_report_upload_opt_in: bool,
    /// When `false` (default), the UI hides adapters whose `tool_license`
    /// doesn't match a standard OSI-approved license. The adapters are
    /// still registered and reachable by id from case.toml — only the
    /// discovery surfaces (wizard / browser tree / command palette /
    /// Tools menu) filter them out.
    ///
    /// Flip to `true` if you need academic-only tools (NAMD, Rosetta,
    /// AlphaFold 3, ChimeraX, VMD, etc.) — see the per-adapter probe
    /// warnings for license terms before commercial use.
    #[serde(default)]
    pub show_non_oss_adapters: bool,
    /// Force the Vina adapter back to subprocess mode (debugging).
    /// When true, overrides any per-case `engine = "native"` choice.
    ///
    /// Wired through to the adapter via
    /// [`valenx_core::set_force_external_vina`], which the app calls
    /// during startup. Provides an escape hatch if the native engine
    /// ever produces wrong results in a real workflow.
    #[serde(default)]
    pub force_external_vina: bool,

    /// User-selected palette variant (Dark / Light / High-Contrast).
    /// Replaces the older [`enum@Theme`] enum's role for variant
    /// selection — that enum still exists for legacy "auto / dark
    /// / light" detection, but the actual painted palette is now
    /// driven by `theme_variant` through [`crate::theme::apply`].
    ///
    /// Defaults to [`ThemeVariant::Dark`] which matches the previous
    /// `Theme::Auto` fallback (dark on platforms without OS-theme
    /// signals). Migration: settings.json files written before
    /// v0.1.1 don't contain this key and serde-default lands them
    /// in `Dark`, which is what `apply_theme` was painting anyway.
    #[serde(default)]
    pub theme_variant: ThemeVariant,

    /// Global font-scale multiplier. Applied at theme-application
    /// time to every entry in `egui::Style::text_styles`. Range
    /// 0.75 – 2.0 (clamped in [`theme::clamp_font_scale`]) — outside
    /// that band, the UI starts losing usability (too small to
    /// click vs popovers off-screen).
    ///
    /// Defaults to `1.0` which leaves egui's baseline sizes
    /// untouched.
    #[serde(default = "default_font_scale")]
    pub font_scale: f32,

    /// Whether the first-launch welcome tour has been shown and
    /// dismissed. Independent of the first-run wizard
    /// (`valenx_first_run::FirstRunDecision`) — the welcome tour is
    /// a 3-step orientation walkthrough that runs on top of the
    /// wizard (which handles adapter detection), not in place of it.
    #[serde(default)]
    pub welcome_tour_completed: bool,

    /// Whether the keyboard-shortcut cheat-sheet overlay should pop
    /// up on the next launch. Toggled by the `?` key + by the Help
    /// menu entry. Persisting it means a user who opened it from
    /// the palette has it remember their preference.
    #[serde(default)]
    pub keyboard_shortcuts_overlay_open: bool,

    /// Most-recently-opened projects, most-recent first. Surfaced on
    /// the welcome landing page so the user can re-open work without
    /// drilling through a folder picker. Capped at
    /// [`Self::MAX_RECENT_PROJECTS`] entries and deduplicated on push;
    /// missing-on-disk entries stay in the list (the landing page
    /// renders them with a "missing" hint rather than silently pruning,
    /// so the user can spot a moved / deleted project).
    ///
    /// Persisted to `settings.json` so the list survives app restarts.
    #[serde(default)]
    pub recent_projects: Vec<std::path::PathBuf>,

    /// Globally suppress every "Open in file browser" action across
    /// the UI. When `true`, all call sites — Settings "Open crashes
    /// folder", Audit "Open audit log location", Run "Open prepared
    /// workdir", Run "Open run workdir", plus the matching command-
    /// palette entries — route through [`crate::file_browser::open_path_or_copy`]
    /// which returns a path-bearing status string instead of spawning
    /// `explorer.exe` / `open` / `xdg-open`. Also suppresses the
    /// `install_crash_reporter` startup-time pre-create of
    /// `<state_dir>/crashes/` — the dir gets created lazily by the
    /// panic hook only when an actual crash is being written.
    ///
    /// **On by default** (since v0.1.0-alpha.2) — a fresh install
    /// never auto-opens File Explorer / Finder / xdg-open for any
    /// reason. Every "Open in file browser" call site instead
    /// displays the path as a status-line message. Users who want
    /// the old "click button → window opens" behaviour untick the
    /// checkbox in Settings → Privacy.
    #[serde(default = "default_disable_file_browser_popups")]
    pub disable_file_browser_popups: bool,

    /// Remember the window's screen position + size across launches so
    /// Valenx reopens where you left it — e.g. on a second monitor. On by
    /// default; untick in Settings → Window to always open at the OS
    /// default spot.
    #[serde(default = "default_remember_window_geometry")]
    pub remember_window_geometry: bool,
    /// Last known window outer position `[x, y]` (physical screen pixels),
    /// persisted while `remember_window_geometry` is on and restored at
    /// startup via the viewport builder. `None` until the window moves.
    #[serde(default)]
    pub window_position: Option<[f32; 2]>,
    /// Last known window inner size `[w, h]` (points), persisted alongside
    /// `window_position`.
    #[serde(default)]
    pub window_size: Option<[f32; 2]>,
    /// Insert a unit starter cube into the viewport whenever a new
    /// project is created (Blender-style). Off by default — CFD / bio
    /// projects don't want stray geometry — so it's an opt-in for
    /// CAD / modelling-first workflows.
    #[serde(default)]
    pub starter_cube_in_new_projects: bool,
}

fn default_disable_file_browser_popups() -> bool {
    true
}

fn default_convergence_target() -> Option<f64> {
    // 1e-4 picks up most reasonably-converged simpleFoam runs without
    // forcing research-grade tolerances. Users tighten or loosen via
    // the Settings dialog.
    Some(1e-4)
}

fn default_font_scale() -> f32 {
    1.0
}

fn default_remember_window_geometry() -> bool {
    true
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            theme: Theme::default(),
            default_shading: ShadingMode::default(),
            residual_scale: ResidualScale::default(),
            convergence_target: default_convergence_target(),
            reprobe_on_close: false,
            crash_report_upload_opt_in: false,
            show_non_oss_adapters: false,
            force_external_vina: false,
            theme_variant: ThemeVariant::default(),
            font_scale: default_font_scale(),
            welcome_tour_completed: false,
            keyboard_shortcuts_overlay_open: false,
            recent_projects: Vec::new(),
            disable_file_browser_popups: default_disable_file_browser_popups(),
            remember_window_geometry: default_remember_window_geometry(),
            window_position: None,
            window_size: None,
            starter_cube_in_new_projects: false,
        }
    }
}

impl Settings {
    /// Cap on the number of recent-project paths surfaced on the
    /// landing page. Keeps the list scannable; the most-recent 8
    /// covers a typical "I'm bouncing between a few cases this week"
    /// workflow without scrolling.
    pub const MAX_RECENT_PROJECTS: usize = 8;

    /// Apply the selected theme variant + font scale to the egui
    /// context. The variant routes through [`theme::apply`] which
    /// builds an egui `Visuals` from the canonical design tokens
    /// (`color::*` in `valenx-design-tokens`). The legacy [`Theme`]
    /// field still influences whether the High-Contrast variant
    /// uses a thicker focus ring — but the actual palette comes
    /// from `theme_variant`.
    pub fn apply_theme(&self, ctx: &egui::Context) {
        theme::apply(ctx, self.theme_variant, self.font_scale);
    }

    /// Push `path` to the front of [`Self::recent_projects`],
    /// deduplicating against any existing entry and trimming the list
    /// to [`Self::MAX_RECENT_PROJECTS`]. The path is canonicalised
    /// best-effort so two opens of the same project via different
    /// relative paths collapse into one entry.
    ///
    /// Returns `true` when the list actually changed (so the caller
    /// can persist) — pushing the already-front entry returns `false`.
    pub fn push_recent_project(&mut self, path: std::path::PathBuf) -> bool {
        // Canonicalise so `./foo/.valenx` and `/abs/foo/.valenx` don't
        // both end up in the list. Falls back to the raw path if the
        // OS reports the file is missing — we still want to remember
        // the user opened it, the landing page just won't be able to
        // re-open it without a fresh pick.
        let canonical = std::fs::canonicalize(&path).unwrap_or(path);
        if self.recent_projects.first() == Some(&canonical) {
            return false;
        }
        self.recent_projects.retain(|p| p != &canonical);
        self.recent_projects.insert(0, canonical);
        self.recent_projects.truncate(Self::MAX_RECENT_PROJECTS);
        true
    }
}

/// Render the Settings window. Returns `true` if any value changed
/// this frame so callers can persist/react.
///
/// The locale catalogue supplies every user-visible label so the
/// Settings panel renders in whichever language the user picked at
/// startup. Pre-i18n call sites that don't have a catalogue handy
/// can use [`show_legacy`] which falls through to the embedded
/// en-US baseline.
/// Action the Settings panel asks the host to perform after the
/// `show` call returns. Lets the panel stay pure-rendering — no
/// filesystem / OS calls in here — while still surfacing
/// user-driven side effects to the caller.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SettingsResponse {
    /// `true` when at least one settings field changed this frame.
    /// Caller persists + reapplies theme based on this flag.
    pub changed: bool,
    /// Set when the user clicked "Open crashes folder…" in the
    /// Privacy section. Caller invokes
    /// [`crate::open_path_in_file_browser`] on the path so the
    /// settings module stays decoupled from the host's
    /// per-platform launcher.
    pub open_crashes_folder: bool,
}

/// Convenience wrapper around [`show_with_response`]; returns the
/// `changed` flag without exposing the full [`SettingsResponse`].
pub fn show(
    ctx: &egui::Context,
    open: &mut bool,
    settings: &mut Settings,
    catalogue: &valenx_i18n::LocaleCatalogue,
) -> bool {
    show_with_response(ctx, open, settings, catalogue).changed
}

/// Like [`show`] but returns the full [`SettingsResponse`] so the
/// caller can react to side-effecting buttons (currently:
/// "Open crashes folder"). Existing call sites that only need
/// the changed-bool can keep using `show`.
pub fn show_with_response(
    ctx: &egui::Context,
    open: &mut bool,
    settings: &mut Settings,
    catalogue: &valenx_i18n::LocaleCatalogue,
) -> SettingsResponse {
    let mut response = SettingsResponse::default();
    let changed = &mut response.changed;

    egui::Window::new(catalogue.lookup("dialog.settings.title"))
        .open(open)
        .collapsible(false)
        .resizable(false)
        .default_width(420.0)
        .show(ctx, |ui| {
            ui.heading(catalogue.lookup("dialog.settings.section.appearance"));
            ui.horizontal(|ui| {
                ui.label(catalogue.lookup("dialog.settings.theme-label"));
                *changed |= ui
                    .radio_value(
                        &mut settings.theme,
                        Theme::Auto,
                        catalogue.lookup("dialog.settings.theme.auto"),
                    )
                    .on_hover_text(
                        "Match the OS preference where supported; falls back to Dark.",
                    )
                    .changed();
                *changed |= ui
                    .radio_value(
                        &mut settings.theme,
                        Theme::Dark,
                        catalogue.lookup("dialog.settings.theme.dark"),
                    )
                    .on_hover_text("Force the dark colour scheme (legacy hint).")
                    .changed();
                *changed |= ui
                    .radio_value(
                        &mut settings.theme,
                        Theme::Light,
                        catalogue.lookup("dialog.settings.theme.light"),
                    )
                    .on_hover_text("Force the light colour scheme (legacy hint).")
                    .changed();
            });
            ui.add_space(4.0);
            // Token-driven palette variant. Replaces the legacy
            // Theme::Auto/Dark/Light split above — that one stays
            // wired so settings.json files from before v0.1.1 keep
            // working, but the actual painted palette comes from
            // `theme_variant` through `theme::apply`.
            ui.horizontal(|ui| {
                ui.label("Palette:");
                for variant in ThemeVariant::ALL {
                    let selected = settings.theme_variant == variant;
                    if ui
                        .radio(selected, variant.label())
                        .on_hover_text(variant.description())
                        .clicked()
                    {
                        settings.theme_variant = variant;
                        *changed = true;
                    }
                }
            });
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                ui.label("Font scale:");
                let mut scale = settings.font_scale;
                if ui
                    .add(
                        egui::Slider::new(&mut scale, 0.75..=2.0)
                            .step_by(0.05)
                            .text("× UI text"),
                    )
                    .on_hover_text(
                        "Multiplier applied to all egui text sizes. 1.0 = default; 1.5 ≈ \
                         WCAG-AAA large text; 2.0 ≈ projector-room readable. Persists.",
                    )
                    .changed()
                {
                    settings.font_scale = theme::clamp_font_scale(scale);
                    *changed = true;
                }
                if ui
                    .small_button("Reset")
                    .on_hover_text("Reset font scale to 1.0 (egui default).")
                    .clicked()
                {
                    settings.font_scale = 1.0;
                    *changed = true;
                }
            });

            ui.add_space(8.0);
            ui.heading(catalogue.lookup("dialog.settings.section.viewport"));
            ui.horizontal(|ui| {
                ui.label(catalogue.lookup("dialog.settings.shading-label"));
                *changed |= ui
                    .radio_value(
                        &mut settings.default_shading,
                        ShadingMode::Shaded,
                        catalogue.lookup("dialog.settings.shading.shaded"),
                    )
                    .changed();
                *changed |= ui
                    .radio_value(
                        &mut settings.default_shading,
                        ShadingMode::Wireframe,
                        catalogue.lookup("dialog.settings.shading.wireframe"),
                    )
                    .changed();
            });

            ui.add_space(8.0);
            ui.heading(catalogue.lookup("dialog.settings.section.residuals"));
            ui.horizontal(|ui| {
                ui.label(catalogue.lookup("dialog.settings.residual-scale-label"));
                *changed |= ui
                    .radio_value(
                        &mut settings.residual_scale,
                        ResidualScale::Log10,
                        catalogue.lookup("dialog.settings.residual-scale.log10"),
                    )
                    .changed();
                *changed |= ui
                    .radio_value(
                        &mut settings.residual_scale,
                        ResidualScale::Linear,
                        catalogue.lookup("dialog.settings.residual-scale.linear"),
                    )
                    .changed();
            });

            ui.add_space(8.0);
            ui.heading("New projects");
            *changed |= ui
                .checkbox(
                    &mut settings.starter_cube_in_new_projects,
                    "Start new projects with a cube (Blender-style)",
                )
                .on_hover_text(
                    "When on, creating a new project drops a unit cube into the \
                     viewport so you have something to model from. Off by default — \
                     analysis projects (CFD, bio) usually want an empty scene.",
                )
                .changed();

            ui.add_space(8.0);
            ui.heading("Window");
            *changed |= ui
                .checkbox(
                    &mut settings.remember_window_geometry,
                    "Remember window position & size (reopen on the same monitor)",
                )
                .on_hover_text(
                    "When on, Valenx reopens where you last left it — including on a \
                     second monitor. Move the window once and it stays there.",
                )
                .changed();
            if ui
                .button("Reset window position")
                .on_hover_text(
                    "Forget the saved position/size so the next launch uses the OS \
                     default — handy if a monitor was unplugged.",
                )
                .clicked()
            {
                settings.window_position = None;
                settings.window_size = None;
                *changed = true;
            }

            ui.add_space(8.0);
            ui.heading(catalogue.lookup("dialog.settings.section.adapters"));
            *changed |= ui
                .checkbox(
                    &mut settings.reprobe_on_close,
                    catalogue.lookup("dialog.settings.reprobe-on-close"),
                )
                .changed();
            // OSS-only filter. Default is OFF (filter on), which
            // hides ~16 adapters wrapping academic-only / non-
            // commercial tools (NAMD, Rosetta, AlphaFold 3,
            // ChimeraX, VMD, NUPACK, mfold, etc.). The adapters
            // stay in the registry — case.toml can still reference
            // them by id — only the discovery surfaces filter them.
            // Plain English here rather than a new locale key so the
            // patch ships without a translation drop.
            *changed |= ui
                .checkbox(
                    &mut settings.show_non_oss_adapters,
                    "Show non-OSS adapters (academic / non-commercial tools)",
                )
                .changed();
            ui.weak(
                "These adapters wrap tools licensed for non-commercial / \
                 academic use only (NAMD, Rosetta, AlphaFold 3, ChimeraX, \
                 VMD, NUPACK, mfold, etc.). The adapters still work — they \
                 just don't appear in the wizard / browser / palette / Tools \
                 menu by default.",
            );
            // Native Vina escape hatch. The adapter normally dispatches
            // to the native engine when a case picks `engine = "native"`;
            // flipping this on calls `valenx_core::set_force_external_vina(true)`
            // at the next launch so the adapter always shells out to
            // the upstream binary instead. Useful for A/B-checking
            // native output against reference Vina.
            *changed |= ui
                .checkbox(
                    &mut settings.force_external_vina,
                    "Force external Vina binary",
                )
                .changed();
            ui.weak(
                "Overrides any per-case `engine = \"native\"` choice and \
                 shells out to the upstream `vina` binary. Takes effect on \
                 next launch.",
            );

            ui.add_space(8.0);
            ui.heading(catalogue.lookup("dialog.settings.section.privacy"));
            *changed |= ui
                .checkbox(
                    &mut settings.crash_report_upload_opt_in,
                    catalogue.lookup("dialog.settings.crash-report-opt-in"),
                )
                .changed();
            ui.weak(catalogue.lookup("dialog.settings.crash-report-explainer"));
            ui.add_space(4.0);
            // "Open crashes folder…" — pure UI signal; the caller
            // resolves the path + invokes the per-OS launcher so
            // this module stays I/O-free.
            let crashes_path = crate::crashes_dir().display().to_string();
            let folder_btn = ui
                .button(catalogue.lookup("dialog.settings.crash-report-open-folder"))
                .on_hover_text(catalogue.format_with(
                    "dialog.settings.crash-report-folder-tooltip",
                    &[("path", &crashes_path)],
                ));
            if folder_btn.clicked() {
                response.open_crashes_folder = true;
            }

            ui.add_space(8.0);
            // Global "no file-browser popups" kill-switch. Off by
            // default — preserves the pre-v0.1.X behaviour where
            // clicking "Open in file browser" actually opens
            // Explorer / Finder / xdg-open. When ON, every such
            // button surfaces the path as a status message instead.
            // Routed through `open_path_or_copy` at every call site.
            // Plain English label rather than a new locale key so the
            // patch ships without a translation drop.
            *changed |= ui
                .checkbox(
                    &mut settings.disable_file_browser_popups,
                    "Never open the system file browser",
                )
                .changed();
            ui.weak(
                "When on, every \"Open in file browser\" button (Settings, \
                 Audit, Run → Open prepared/run workdir, command palette) \
                 shows the path as a status message instead of launching \
                 Explorer / Finder / xdg-open. Also stops pre-creating the \
                 crashes folder on startup — it appears only when an actual \
                 crash report is written.",
            );

            ui.add_space(12.0);
            ui.separator();
            ui.horizontal(|ui| {
                if ui
                    .button(catalogue.lookup("dialog.settings.reset-defaults"))
                    .clicked()
                {
                    *settings = Settings::default();
                    *changed = true;
                }
                ui.weak(catalogue.lookup("dialog.settings.persistence-note"));
            });
        });

    response
}

/// Backwards-compatible shim that loads the embedded en-US
/// catalogue and dispatches to [`show`]. Lets pre-i18n call sites
/// keep their `show(ctx, open, settings)` signature while the
/// rest of the workspace migrates over.
pub fn show_legacy(ctx: &egui::Context, open: &mut bool, settings: &mut Settings) -> bool {
    let catalogue = valenx_i18n::embedded_en_us();
    show(ctx, open, settings, &catalogue)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_sane() {
        let s = Settings::default();
        assert_eq!(s.theme, Theme::Auto);
        assert_eq!(s.residual_scale, ResidualScale::Log10);
        assert_eq!(s.default_shading, ShadingMode::Shaded);
        assert!(!s.reprobe_on_close);
        // Privacy: opt-in upload is OFF by default. The reporter
        // still writes reports to disk; this flag only governs the
        // network egress on the next launch's "submit?" prompt.
        assert!(!s.crash_report_upload_opt_in);
        // OSS-only filter is ON by default — non-OSS adapters are
        // hidden from the wizard / browser / palette / Tools menu.
        // Users opt back in via the Settings dialog when they need
        // academic-only tools (NAMD, Rosetta, AlphaFold 3, etc.).
        assert!(!s.show_non_oss_adapters);
        // Polish-pass defaults: the token-driven palette starts in
        // Dark (matching the legacy auto-on-no-OS-signal behaviour),
        // the font scale at 1.0 (egui defaults), and neither the
        // welcome tour nor the shortcuts overlay open on launch
        // (the welcome tour is gated by a separate first-run check
        // on the next launch after install).
        assert_eq!(s.theme_variant, ThemeVariant::Dark);
        assert_eq!(s.font_scale, 1.0);
        assert!(!s.welcome_tour_completed);
        assert!(!s.keyboard_shortcuts_overlay_open);
        // The file-browser-popup kill-switch defaults ON since
        // v0.1.0-alpha.2 — a fresh install never auto-opens File
        // Explorer / Finder / xdg-open. Users who want the old
        // "click button → window opens" behaviour untick the
        // checkbox in Settings → Privacy.
        assert!(s.disable_file_browser_popups);
    }

    #[test]
    fn theme_variant_round_trips_through_serde() {
        let s = Settings {
            theme_variant: ThemeVariant::HighContrast,
            font_scale: 1.5,
            welcome_tour_completed: true,
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let back: Settings = serde_json::from_str(&json).unwrap();
        assert_eq!(back.theme_variant, ThemeVariant::HighContrast);
        assert_eq!(back.font_scale, 1.5);
        assert!(back.welcome_tour_completed);
    }

    #[test]
    fn new_settings_fields_default_when_missing_from_json() {
        // Older settings.json (written before the polish pass) lacks
        // every new field — serde-defaults must produce the same
        // values as a fresh `Settings::default()`.
        let json = r#"{"theme":"auto","default_shading":"shaded","residual_scale":"log10","reprobe_on_close":false}"#;
        let s: Settings = serde_json::from_str(json).unwrap();
        assert_eq!(s.theme_variant, ThemeVariant::Dark);
        assert_eq!(s.font_scale, 1.0);
        assert!(!s.welcome_tour_completed);
        assert!(!s.keyboard_shortcuts_overlay_open);
    }

    #[test]
    fn show_non_oss_adapters_round_trips_through_serde() {
        let s = Settings {
            show_non_oss_adapters: true,
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let back: Settings = serde_json::from_str(&json).unwrap();
        assert!(back.show_non_oss_adapters);
    }

    #[test]
    fn show_non_oss_adapters_defaults_off_when_field_missing_in_json() {
        // Old settings.json files (saved before this field existed)
        // must continue to load with the OSS-only filter ON.
        let json = r#"{"theme":"auto","default_shading":"shaded","residual_scale":"log10","reprobe_on_close":false}"#;
        let s: Settings = serde_json::from_str(json).unwrap();
        assert!(!s.show_non_oss_adapters);
    }

    #[test]
    fn crash_report_opt_in_round_trips_through_serde() {
        let s = Settings {
            crash_report_upload_opt_in: true,
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let back: Settings = serde_json::from_str(&json).unwrap();
        assert!(back.crash_report_upload_opt_in);
    }

    #[test]
    fn crash_report_opt_in_defaults_off_when_field_missing_in_json() {
        // serde-default fires on deserialisation when the JSON
        // doesn't have the key. Mirrors what an old settings.json
        // file looks like after upgrading.
        let json = r#"{"theme":"auto","default_shading":"shaded","residual_scale":"log10","reprobe_on_close":false}"#;
        let s: Settings = serde_json::from_str(json).unwrap();
        assert!(!s.crash_report_upload_opt_in);
    }

    #[test]
    fn convergence_target_defaults_to_engineering_value() {
        // Default 1e-4 catches most reasonably-converged simpleFoam
        // runs without forcing research-grade tolerances. Bumping
        // this default should be a deliberate edit visible in
        // review.
        let s = Settings::default();
        assert_eq!(s.convergence_target, Some(1e-4));
    }

    #[test]
    fn convergence_target_round_trips_through_serde() {
        let mut s = Settings {
            convergence_target: Some(1e-6),
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let parsed: Settings = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.convergence_target, Some(1e-6));

        // None round-trips too (user disabled the overlay).
        s.convergence_target = None;
        let json = serde_json::to_string(&s).unwrap();
        let parsed: Settings = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.convergence_target, None);
    }

    #[test]
    fn recent_projects_defaults_empty() {
        // Brand-new install has no project history. The landing page
        // renders the "(no recent projects yet — start by creating
        // one)" hint when this is empty.
        let s = Settings::default();
        assert!(s.recent_projects.is_empty());
    }

    #[test]
    fn push_recent_project_dedups_and_moves_to_front() {
        // Pushing the same path twice must NOT produce two entries; the
        // existing entry moves to the front instead. Use a path the
        // OS won't canonicalise (a fictitious deep path) so the test
        // doesn't depend on the temp-dir being real.
        let mut s = Settings::default();
        let a = std::path::PathBuf::from("/no/such/dir/a.valenx");
        let b = std::path::PathBuf::from("/no/such/dir/b.valenx");
        assert!(s.push_recent_project(a.clone()));
        assert!(s.push_recent_project(b.clone()));
        assert!(s.push_recent_project(a.clone()));
        assert_eq!(s.recent_projects.len(), 2);
        assert_eq!(s.recent_projects[0], a);
        assert_eq!(s.recent_projects[1], b);
    }

    #[test]
    fn push_recent_project_caps_at_max() {
        // Pushing past the cap evicts the oldest entries.
        let mut s = Settings::default();
        for i in 0..(Settings::MAX_RECENT_PROJECTS + 4) {
            s.push_recent_project(std::path::PathBuf::from(format!(
                "/no/such/dir/p-{i}.valenx"
            )));
        }
        assert_eq!(s.recent_projects.len(), Settings::MAX_RECENT_PROJECTS);
        // The first entry is the LAST pushed (most-recent-first order).
        assert!(s.recent_projects[0]
            .to_string_lossy()
            .ends_with(&format!(
                "p-{}.valenx",
                Settings::MAX_RECENT_PROJECTS + 3
            )));
    }

    #[test]
    fn push_recent_project_returns_false_when_already_front() {
        // Pushing the already-front entry is a no-op — the dedup
        // branch returns `false` so the caller skips disk persistence.
        let mut s = Settings::default();
        let a = std::path::PathBuf::from("/no/such/dir/a.valenx");
        assert!(s.push_recent_project(a.clone()));
        assert!(!s.push_recent_project(a));
    }

    #[test]
    fn recent_projects_round_trips_through_serde() {
        let s = Settings {
            recent_projects: vec![
                std::path::PathBuf::from("/no/such/dir/a.valenx"),
                std::path::PathBuf::from("/no/such/dir/b.valenx"),
            ],
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let back: Settings = serde_json::from_str(&json).unwrap();
        assert_eq!(back.recent_projects.len(), 2);
    }

    #[test]
    fn disable_file_browser_popups_round_trips_through_serde() {
        let s = Settings {
            disable_file_browser_popups: true,
            ..Settings::default()
        };
        let json = serde_json::to_string(&s).unwrap();
        let back: Settings = serde_json::from_str(&json).unwrap();
        assert!(back.disable_file_browser_popups);
    }

    #[test]
    fn disable_file_browser_popups_defaults_on_when_field_missing_in_json() {
        // Older settings.json files (saved before this field existed)
        // load with the kill-switch ON since v0.1.0-alpha.2 — a fresh
        // install never auto-opens File Explorer / Finder / xdg-open.
        // Users who want the old behaviour untick the checkbox in
        // Settings → Privacy.
        let json = r#"{"theme":"auto","default_shading":"shaded","residual_scale":"log10","reprobe_on_close":false}"#;
        let s: Settings = serde_json::from_str(json).unwrap();
        assert!(s.disable_file_browser_popups);
    }

    #[test]
    fn push_recent_project_canonicalises_real_dirs() {
        // The dedup branch only collapses two paths if `canonicalize`
        // succeeds and returns the same absolute path for both. With
        // fictitious paths the previous tests can't reach that branch
        // (canonicalize errors out and the unwrap_or falls back to
        // the raw input). Exercise the real branch by pushing the
        // same physical directory via two different representations
        // (an absolute form and a `./`-prefixed relative form), and
        // assert the list collapses to one entry.
        let tmp = tempfile::tempdir().expect("create tempdir");
        let abs = std::fs::canonicalize(tmp.path()).expect("canonicalize tempdir");

        // Push absolute path first.
        let mut s = Settings::default();
        assert!(s.push_recent_project(abs.clone()));
        assert_eq!(s.recent_projects.len(), 1);

        // Push the same directory via a `./`-prefixed relative path
        // resolved against the tempdir's parent. We chdir into the
        // parent so `./<leaf>` references the same physical directory
        // as `abs`. The canonicalisation inside `push_recent_project`
        // must collapse the two representations.
        let parent = abs.parent().expect("tempdir has a parent");
        let leaf = abs.file_name().expect("tempdir has a leaf");
        let prev_cwd = std::env::current_dir().expect("read cwd");
        std::env::set_current_dir(parent).expect("chdir to parent");
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let relative = std::path::PathBuf::from(".").join(leaf);
            // The second push should be a no-op (already at front
            // after canonicalisation collapses to the same `abs`).
            // The first invariant — list length stays at 1 — is the
            // canonical dedup proof.
            let _changed = s.push_recent_project(relative);
            assert_eq!(
                s.recent_projects.len(),
                1,
                "canonicalisation should collapse `./leaf` and abs into one entry"
            );
            assert_eq!(s.recent_projects[0], abs);
        }));
        std::env::set_current_dir(prev_cwd).expect("restore cwd");
        result.expect("inner assertions");
    }
}
