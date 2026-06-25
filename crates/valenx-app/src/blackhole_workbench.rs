//! The right-side **Black-Hole / Relativity Workbench** panel — a native
//! front-end over the in-house `valenx-relativity` general-relativity engine.
//!
//! Mirrors the other workbenches (`fem_workbench`, `aero_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! `crate::ValenxApp::show_blackhole_workbench`, toggled from the View menu.
//!
//! v1 surfaces the (already complete + closed-form-validated) observables and
//! thermodynamics of the **Kerr–Newman family** — pick a spacetime and set
//! `(M, a, Q)`, then read horizons / photon sphere / ISCO / shadow radius /
//! gravitational redshift / Hawking temperature / entropy / horizon angular
//! velocity — plus a backward-ray-traced **shadow / lensing image** rendered by
//! `valenx_relativity::shadow::render_shadow` and shown as an egui texture (the
//! same `ColorImage::from_rgb → load_texture` path as `render_workbench.rs`).
//!
//! Honesty: a super-extremal hole (`a² + Q² > M²`, a naked singularity) or a
//! quantity with no implemented closed form (e.g. an ISCO with both spin and
//! charge) makes the engine return a [`RelativityError`]; the workbench surfaces
//! that error verbatim rather than inventing a number. No new physics — wiring
//! and presentation over the validated engine only.

use std::fmt::Write as _;

use eframe::egui;

use valenx_relativity::shadow::{render_shadow, ShadowImage};
use valenx_relativity::spacetimes::KerrNewman;
use valenx_relativity::{
    gravitational_redshift, horizons, isco, photon_sphere, shadow_radius, thermodynamics,
    OrbitSense, RelativityError,
};

use crate::ValenxApp;

/// Which named member of the Kerr–Newman family the user is exploring.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SpacetimeKind {
    /// Non-rotating, uncharged (`a = Q = 0`).
    #[default]
    Schwarzschild,
    /// Rotating, uncharged (`Q = 0`).
    Kerr,
    /// Charged, non-rotating (`a = 0`).
    ReissnerNordstrom,
    /// The general rotating + charged hole.
    KerrNewman,
}

impl SpacetimeKind {
    /// Human-readable label for the picker.
    fn label(self) -> &'static str {
        match self {
            SpacetimeKind::Schwarzschild => "Schwarzschild (M)",
            SpacetimeKind::Kerr => "Kerr (M, a)",
            SpacetimeKind::ReissnerNordstrom => "Reissner–Nordström (M, Q)",
            SpacetimeKind::KerrNewman => "Kerr–Newman (M, a, Q)",
        }
    }

    /// Whether this member has a spin degree of freedom (enables the `a` field).
    fn has_spin(self) -> bool {
        matches!(self, SpacetimeKind::Kerr | SpacetimeKind::KerrNewman)
    }

    /// Whether this member has a charge degree of freedom (enables the `Q` field).
    fn has_charge(self) -> bool {
        matches!(
            self,
            SpacetimeKind::ReissnerNordstrom | SpacetimeKind::KerrNewman
        )
    }
}

/// Persistent state for the Black-Hole / Relativity workbench.
pub struct BlackHoleWorkbenchState {
    /// Which member of the Kerr–Newman family is selected.
    pub kind: SpacetimeKind,
    /// Mass `M` (geometrized units; the natural scale is `M = 1`).
    pub mass: f64,
    /// Spin parameter `a = J/M`.
    pub spin: f64,
    /// Electric charge `Q`.
    pub charge: f64,
    /// Observer radius (units of `M`) for the shadow render.
    pub r_obs: f64,
    /// Observer polar angle (degrees) for the shadow render.
    pub theta_obs_deg: f64,
    /// Shadow image resolution (pixels, square).
    pub img_size: usize,
    /// Formatted observables readout (or a `⚠ <error>` string).
    pub readout: String,
    /// The rendered shadow texture, if any.
    pub texture: Option<egui::TextureHandle>,
}

impl Default for BlackHoleWorkbenchState {
    fn default() -> Self {
        Self {
            kind: SpacetimeKind::Schwarzschild,
            mass: 1.0,
            spin: 0.0,
            charge: 0.0,
            r_obs: 50.0,
            theta_obs_deg: 80.0,
            img_size: 160,
            readout: String::new(),
            texture: None,
        }
    }
}

/// Parse a spacetime-family name (for the agent `SetControl` bridge) into a
/// [`SpacetimeKind`]. Case-insensitive; accepts short names + the picker words.
/// Fail-loud on an unknown name.
fn parse_spacetime_kind(s: &str) -> Result<SpacetimeKind, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "schwarzschild" => Ok(SpacetimeKind::Schwarzschild),
        "kerr" => Ok(SpacetimeKind::Kerr),
        "reissner-nordstrom" | "reissner-nordström" | "reissnernordstrom" | "rn" => {
            Ok(SpacetimeKind::ReissnerNordstrom)
        }
        "kerr-newman" | "kerrnewman" | "kn" => Ok(SpacetimeKind::KerrNewman),
        other => Err(format!(
            "unknown spacetime '{other}' (expected 'schwarzschild', 'kerr', \
             'reissner-nordstrom', or 'kerr-newman')"
        )),
    }
}

impl BlackHoleWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). `spacetime` is the family
    /// picker (set by option name); the rest are numeric.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "spacetime",
            "mass M",
            "spin a",
            "charge Q",
            "observer r (M)",
            "observer \u{03B8} (\u{00B0})",
            "image size (px)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Numeric fields read [`AgentValue::as_f64`] /
    /// [`AgentValue::as_i64`]; `spacetime` reads [`AgentValue::as_str`] and
    /// matches the family names. Fail-loud: an unknown caption, wrong type, or
    /// negative image size returns `Err` — never a panic, no field written on
    /// error. (Setting `spin`/`charge` on a member that ignores them is allowed;
    /// [`hole`](Self::hole) zeroes the unused degree of freedom at solve time.)
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "spacetime" => self.kind = parse_spacetime_kind(value.as_str()?)?,
            "mass M" => self.mass = value.as_f64()?,
            "spin a" => self.spin = value.as_f64()?,
            "charge Q" => self.charge = value.as_f64()?,
            "observer r (M)" => self.r_obs = value.as_f64()?,
            "observer \u{03B8} (\u{00B0})" => self.theta_obs_deg = value.as_f64()?,
            "image size (px)" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("image size must be >= 0, got {n}"));
                }
                self.img_size = n as usize;
            }
            other => return Err(format!("unknown Black Hole control: {other:?}")),
        }
        Ok(())
    }

    /// Build the parameter struct, zeroing spin/charge for the simpler members
    /// so the UI can't produce a contradictory hole.
    pub fn hole(&self) -> KerrNewman {
        match self.kind {
            SpacetimeKind::Schwarzschild => KerrNewman {
                mass: self.mass,
                spin: 0.0,
                charge: 0.0,
            },
            SpacetimeKind::Kerr => KerrNewman {
                mass: self.mass,
                spin: self.spin,
                charge: 0.0,
            },
            SpacetimeKind::ReissnerNordstrom => KerrNewman {
                mass: self.mass,
                spin: 0.0,
                charge: self.charge,
            },
            SpacetimeKind::KerrNewman => KerrNewman {
                mass: self.mass,
                spin: self.spin,
                charge: self.charge,
            },
        }
    }
}

/// Radius (in `M`) at which the gravitational-redshift "observer at infinity" is
/// evaluated. The engine takes a finite `r_obs`; a far radius makes `−g_tt`
/// effectively unity (`1 − 2M/r → 1`) while staying a real static observer.
const REDSHIFT_OBSERVER_RADIUS: f64 = 1.0e6;

/// Format the full observables readout, or return the engine's error verbatim.
///
/// A super-extremal hole (`a² + Q² > M²`, a naked singularity), a non-positive
/// mass, or a quantity with no implemented closed form (e.g. ISCO/photon-sphere
/// with both spin and charge, or a closed-form shadow radius for a *rotating*
/// hole) surfaces honestly as a [`RelativityError`] — the workbench never fakes
/// a number.
///
/// # Errors
/// Propagates any [`RelativityError`] from the engine observables/thermo calls.
pub fn compute_observables(s: &BlackHoleWorkbenchState) -> Result<String, RelativityError> {
    let bh = s.hole();
    // Validate up front (also what every observable does internally): this makes
    // a super-extremal hole fail before we format a partial readout.
    let h = horizons(&bh)?;
    let thermo = thermodynamics(&bh)?;

    let mut out = String::new();
    let _ = writeln!(out, "spacetime         : {}", s.kind.label());
    let _ = writeln!(out, "mass M            : {:.3}", bh.mass);
    let _ = writeln!(out, "spin a            : {:.3}", bh.spin);
    let _ = writeln!(out, "charge Q          : {:.3}", bh.charge);
    let _ = writeln!(out, "event horizon r+  : {:.3} M", h.outer);
    let _ = writeln!(out, "inner horizon r-  : {:.3} M", h.inner);

    // Photon sphere / ISCO / shadow radius / redshift each have parameter ranges
    // with no closed form; show the value where available, else the engine's
    // reason (Unsupported …) on that line rather than aborting the whole card.
    match photon_sphere(&bh, OrbitSense::Prograde) {
        Ok(r) => {
            let _ = writeln!(out, "photon sphere     : {r:.3} M (prograde)");
        }
        Err(e) => {
            let _ = writeln!(out, "photon sphere     : — ({e})");
        }
    }
    match isco(&bh, OrbitSense::Prograde) {
        Ok(r) => {
            let _ = writeln!(out, "ISCO (prograde)   : {r:.3} M");
            // Representative gravitational redshift: light from the ISCO seen by
            // a distant static observer.
            match gravitational_redshift(&bh, r, REDSHIFT_OBSERVER_RADIUS) {
                Ok(z) => {
                    let _ = writeln!(out, "redshift 1+z @ISCO: {z:.4}");
                }
                Err(e) => {
                    let _ = writeln!(out, "redshift 1+z @ISCO: — ({e})");
                }
            }
        }
        Err(e) => {
            let _ = writeln!(out, "ISCO (prograde)   : — ({e})");
        }
    }
    match shadow_radius(&bh) {
        Ok(r) => {
            let _ = writeln!(out, "shadow radius     : {r:.3} M");
        }
        Err(e) => {
            let _ = writeln!(out, "shadow radius     : — ({e})");
        }
    }

    let _ = writeln!(
        out,
        "Hawking T_H       : {:.4e}",
        thermo.hawking_temperature
    );
    let _ = writeln!(out, "entropy S         : {:.4e}", thermo.entropy);
    let _ = writeln!(
        out,
        "horizon Ω_H       : {:.4}",
        thermo.horizon_angular_velocity
    );
    Ok(out)
}

/// Map a shadow mask to an RGB8 buffer (row-major, `w*h*3` bytes): the captured
/// region is black; escaping pixels get a simple radial "sky" tint so the photon
/// ring reads as the boundary. Cosmetic only — richer shading (an accretion
/// disk) is a later phase.
pub fn shadow_to_rgb(img: &ShadowImage) -> Vec<u8> {
    let (w, h) = (img.width, img.height);
    let mut rgb = vec![0u8; w * h * 3];
    let (cx, cy) = (w as f64 / 2.0, h as f64 / 2.0);
    let norm = cx.max(cy).max(1.0);
    for row in 0..h {
        for col in 0..w {
            let i = (row * w + col) * 3;
            if img.is_shadow(col, row) {
                continue; // captured by the hole — black
            }
            // Sky tint by distance from image centre (purely cosmetic).
            let d = (((col as f64 - cx).powi(2) + (row as f64 - cy).powi(2)).sqrt() / norm)
                .clamp(0.0, 1.0);
            rgb[i] = (40.0 + 120.0 * d) as u8;
            rgb[i + 1] = (60.0 + 120.0 * d) as u8;
            rgb[i + 2] = (120.0 + 135.0 * d) as u8;
        }
    }
    rgb
}

/// Image-plane half-extent (in `M`) for the shadow render: a few × the
/// Schwarzschild shadow radius (`3√3 M ≈ 5.2 M`) so the photon ring sits well
/// inside the frame for the spins/charges the UI allows.
const SHADOW_HALF_EXTENT: f64 = 12.0;

/// Render the shadow for the current state into an egui texture.
///
/// # Errors
/// Propagates [`RelativityError`] from the ray-tracer (e.g. an observer too
/// close to the hole, or on the polar axis).
pub fn render_shadow_texture(
    s: &BlackHoleWorkbenchState,
    ctx: &egui::Context,
) -> Result<egui::TextureHandle, RelativityError> {
    let bh = s.hole();
    let n = s.img_size.clamp(16, 256);
    let img: ShadowImage = render_shadow(
        &bh,
        s.r_obs,
        s.theta_obs_deg.to_radians(),
        SHADOW_HALF_EXTENT,
        n,
        n,
    )?;
    let rgb = shadow_to_rgb(&img);
    let color = egui::ColorImage::from_rgb([img.width, img.height], &rgb);
    Ok(ctx.load_texture("blackhole_shadow", color, egui::TextureOptions::LINEAR))
}

/// Draw the Black-Hole / Relativity workbench (a no-op unless toggled on via
/// View → Black Hole). Mirrors [`crate::fem_workbench::draw_fem_workbench`].
pub fn draw_blackhole_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_blackhole_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_blackhole_workbench",
        "Black Hole",
        blackhole_workbench_body,
    );
    if close {
        app.show_blackhole_workbench = false;
    }
}

/// The workbench body: spacetime picker + `(M, a, Q)` + observer controls, a
/// "Compute observables" button, a "Render shadow" button, and the readout +
/// shadow image. Split out so [`workbench_shell`](crate::workbench_chrome::workbench_shell)
/// (and any future dockable layout) can call it with the same `(app, ui)` shape.
fn blackhole_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("general relativity · valenx-relativity")
            .weak()
            .small(),
    );
    ui.separator();

    let mut do_compute = false;
    let mut do_render = false;
    {
        let s = &mut app.blackhole;

        egui::ComboBox::from_label("spacetime")
            .selected_text(s.kind.label())
            .show_ui(ui, |ui| {
                for kind in [
                    SpacetimeKind::Schwarzschild,
                    SpacetimeKind::Kerr,
                    SpacetimeKind::ReissnerNordstrom,
                    SpacetimeKind::KerrNewman,
                ] {
                    ui.selectable_value(&mut s.kind, kind, kind.label());
                }
            });

        // Associate each numeric `DragValue` with its grid caption via
        // `labelled_by`, so the spin button carries the caption as its
        // accessibility / UI-Automation Name (egui clears a DragValue's own
        // Name, leaving it anonymous to a screen reader / AI driver otherwise);
        // the hover text mirrors the caption for a mouse user.
        egui::Grid::new("blackhole_params")
            .num_columns(2)
            .show(ui, |ui| {
                let m = ui.label("mass M");
                ui.add(
                    egui::DragValue::new(&mut s.mass)
                        .speed(0.05)
                        .range(0.01..=1.0e3),
                )
                .labelled_by(m.id)
                .on_hover_text("Mass M (geometrized units)");
                ui.end_row();

                let a = ui.label("spin a");
                ui.add_enabled(
                    s.kind.has_spin(),
                    egui::DragValue::new(&mut s.spin)
                        .speed(0.01)
                        .range(0.0..=1.0e3),
                )
                .labelled_by(a.id)
                .on_hover_text("Spin parameter a = J/M");
                ui.end_row();

                let q = ui.label("charge Q");
                ui.add_enabled(
                    s.kind.has_charge(),
                    egui::DragValue::new(&mut s.charge)
                        .speed(0.01)
                        .range(0.0..=1.0e3),
                )
                .labelled_by(q.id)
                .on_hover_text("Electric charge Q");
                ui.end_row();
            });

        if ui
            .button(egui::RichText::new("Compute observables").strong())
            .clicked()
        {
            do_compute = true;
        }

        ui.add_space(6.0);
        ui.label(egui::RichText::new("Shadow / lensing image").strong());
        // Illustrative presets (geometrized units; spins are representative,
        // shown as a*=a/M for the famous EHT targets — not fitted parameters).
        ui.horizontal(|ui| {
            if ui
                .button("Sgr A*")
                .on_hover_text("preset: Kerr, M = 1, a ≈ 0.50 (illustrative)")
                .clicked()
            {
                s.kind = SpacetimeKind::Kerr;
                s.mass = 1.0;
                s.spin = 0.50;
                s.charge = 0.0;
            }
            if ui
                .button("M87*")
                .on_hover_text("preset: Kerr, M = 1, a ≈ 0.90 (illustrative)")
                .clicked()
            {
                s.kind = SpacetimeKind::Kerr;
                s.mass = 1.0;
                s.spin = 0.90;
                s.charge = 0.0;
            }
        });

        egui::Grid::new("blackhole_observer")
            .num_columns(2)
            .show(ui, |ui| {
                let r = ui.label("observer r (M)");
                ui.add(
                    egui::DragValue::new(&mut s.r_obs)
                        .speed(1.0)
                        .range(10.0..=1.0e4),
                )
                .labelled_by(r.id)
                .on_hover_text("Observer radius r (units of M)");
                ui.end_row();
                let theta = ui.label("observer θ (°)");
                ui.add(
                    egui::DragValue::new(&mut s.theta_obs_deg)
                        .speed(1.0)
                        .range(1.0..=179.0),
                )
                .labelled_by(theta.id)
                .on_hover_text("Observer polar angle θ (degrees)");
                ui.end_row();
                let px = ui.label("image size (px)");
                ui.add(
                    egui::DragValue::new(&mut s.img_size)
                        .speed(4.0)
                        .range(16..=256),
                )
                .labelled_by(px.id)
                .on_hover_text("Shadow image resolution (pixels, square)");
                ui.end_row();
            });

        if ui.button("Render shadow").clicked() {
            do_render = true;
        }
    }

    // Run the (synchronous) compute/render outside the field borrow above so we
    // can also borrow `ui.ctx()`. The shadow ray-trace is capped at 256² so a
    // click can't wedge the frame; a background-thread render for larger sizes
    // is a later pass (cf. render_workbench's BackgroundJob).
    if do_compute {
        app.blackhole.readout = match compute_observables(&app.blackhole) {
            Ok(t) => t,
            Err(e) => format!("⚠ {e}"),
        };
    }
    if do_render {
        match render_shadow_texture(&app.blackhole, ui.ctx()) {
            Ok(tex) => app.blackhole.texture = Some(tex),
            Err(e) => {
                app.blackhole.readout = format!("⚠ {e}");
                app.blackhole.texture = None;
            }
        }
    }

    let s = &app.blackhole;
    if !s.readout.is_empty() {
        ui.add_space(6.0);
        ui.add(
            egui::TextEdit::multiline(&mut s.readout.as_str())
                .font(egui::TextStyle::Monospace)
                .desired_width(f32::INFINITY)
                .interactive(false),
        );
    }
    if let Some(tex) = &s.texture {
        ui.add_space(4.0);
        ui.add(
            egui::Image::new(egui::load::SizedTexture::new(tex.id(), tex.size_vec2()))
                .max_width(ui.available_width()),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        use crate::agent_commands::AgentValue;
        let mut s = BlackHoleWorkbenchState::default();
        // A representative float set lands in state.
        s.agent_set("mass M", &AgentValue::Float(2.5)).unwrap();
        assert_eq!(s.mass, 2.5);
        // The spacetime enum is set by option name.
        s.agent_set("spacetime", &AgentValue::Str("kerr".into()))
            .unwrap();
        assert_eq!(s.kind, SpacetimeKind::Kerr);
        // The non-ASCII observer-angle caption resolves.
        s.agent_set("observer \u{03B8} (\u{00B0})", &AgentValue::Float(45.0))
            .unwrap();
        assert_eq!(s.theta_obs_deg, 45.0);
        // Unknown caption -> Err.
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into the float mass) -> Err, field untouched.
        assert!(s
            .agent_set("mass M", &AgentValue::Str("heavy".into()))
            .is_err());
        assert_eq!(s.mass, 2.5, "rejected set leaves field untouched");
    }

    #[test]
    fn schwarzschild_observables_match_closed_form() {
        let s = BlackHoleWorkbenchState {
            mass: 1.0,
            ..Default::default()
        };
        let text = compute_observables(&s).expect("valid hole");
        // Event horizon r+ = 2M, photon sphere 3M, ISCO 6M for M = 1.
        assert!(text.contains("2.000"), "horizon r+=2M:\n{text}");
        assert!(text.contains("3.000"), "photon sphere 3M:\n{text}");
        assert!(text.contains("6.000"), "ISCO 6M:\n{text}");
    }

    #[test]
    fn super_extremal_is_rejected_not_faked() {
        // a² + Q² > M² is a naked singularity — must error, not invent numbers.
        let s = BlackHoleWorkbenchState {
            kind: SpacetimeKind::Kerr,
            mass: 1.0,
            spin: 2.0,
            ..Default::default()
        };
        assert!(compute_observables(&s).is_err());
    }

    #[test]
    fn shadow_mask_maps_to_black_inside_and_nonblack_outside() {
        let img = ShadowImage {
            width: 2,
            height: 1,
            half_extent: 8.0,
            shadow: vec![true, false], // pixel 0 captured, pixel 1 escapes
        };
        let rgb = shadow_to_rgb(&img);
        assert_eq!(rgb.len(), 2 * 1 * 3);
        assert_eq!(&rgb[0..3], &[0, 0, 0], "shadow pixel is black");
        assert!(rgb[3..6].iter().any(|&c| c > 0), "sky pixel is not black");
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// Draw the panel once in a headless egui context **with accesskit enabled**
    /// and return the emitted accessibility tree nodes — the same tree a screen
    /// reader / AI UI-Automation driver consumes. `accesskit` is re-exported by
    /// egui, so this needs no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_blackhole_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_blackhole_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_blackhole_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_blackhole_workbench = true;
        // A populated readout must also draw without panic.
        app.blackhole.readout = compute_observables(&app.blackhole).unwrap();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The (M, a, Q) and observer DragValues are SpinButtons; each must be
        // `labelled_by` its grid caption (egui clears a DragValue's own Name) —
        // including the spin/charge fields, which are present in the tree even
        // when disabled for the Schwarzschild default.
        let mut app = ValenxApp::default();
        app.show_blackhole_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // mass, spin, charge, observer r, observer θ, image size = 6.
        assert!(
            spin_buttons.len() >= 6,
            "expected the black-hole numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every black-hole DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["mass M", "spin a", "charge Q", "observer r (M)"] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The spacetime ComboBox and the action buttons stay named + invokable.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::ComboBox),
            "the spacetime picker is a ComboBox node"
        );
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Compute observables"))),
            "the Compute observables button is a named, invokable node"
        );
    }
}
