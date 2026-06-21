//! The right-side **Solenoid Coil Workbench** panel — native long-solenoid
//! inductor analysis over `valenx-coil`.
//!
//! Mirrors the Heat Transfer / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_coil_workbench`,
//! toggled from the View menu. The form sets a uniformly wound solenoid
//! (turns `N`, radius, length, current, core material) and "Analyze"
//! reports the axial field, self-inductance, ampere-turns, stored energy
//! and the inductive reactance from the closed-form long-solenoid model;
//! "Show 3-D coil" loads a representative winding cylinder (with a thin
//! inner core rod) into the central viewport.

use std::f64::consts::{PI, TAU};
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_coil::{energy_joules, reactance_ohms, Solenoid, VACUUM_PERMEABILITY};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The core material wound through the solenoid, selecting the relative
/// permeability `mu_r` used in `L = mu0 * mu_r * N^2 * A / l`. Representative
/// textbook figures; real cores are frequency- and drive-dependent.
#[derive(Debug, Clone, Copy, PartialEq)]
enum CoreKind {
    /// Air / vacuum core, `mu_r = 1`.
    Air,
    /// Soft-iron core, representative `mu_r ~ 200`.
    Iron,
    /// Ferrite core, representative `mu_r ~ 2000`.
    Ferrite,
}

impl CoreKind {
    /// The relative permeability `mu_r` (dimensionless, `>= 1`) this core
    /// contributes to the inductance.
    fn relative_permeability(self) -> f64 {
        match self {
            CoreKind::Air => 1.0,
            CoreKind::Iron => 200.0,
            CoreKind::Ferrite => 2000.0,
        }
    }

    /// A short human label for the readout.
    fn label(self) -> &'static str {
        match self {
            CoreKind::Air => "air",
            CoreKind::Iron => "iron",
            CoreKind::Ferrite => "ferrite",
        }
    }
}

/// Persistent form + result state for the Solenoid Coil Workbench.
pub struct CoilWorkbenchState {
    /// Number of turns `N` (dimensionless).
    turns: f64,
    /// Winding radius `r` (m); the cross-sectional area is `pi * r^2`.
    radius_m: f64,
    /// Axial winding length `l` (m).
    length_m: f64,
    /// Steady coil current `I` (A).
    current_a: f64,
    /// Drive frequency `f` (Hz) for the inductive reactance.
    frequency_hz: f64,
    /// Series circuit resistance `R` (ohm) for the RL time constant
    /// `tau = L / R`.
    series_resistance_ohms: f64,
    /// Core material, selecting the relative permeability `mu_r`.
    core: CoreKind,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D coil solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for CoilWorkbenchState {
    fn default() -> Self {
        // A 500-turn air-cored solenoid, 1 cm winding radius, 10 cm long,
        // carrying 2 A at 1 kHz with 4 ohm of series resistance:
        // B ~ 12.6 mT, L ~ 0.99 mH, X_L ~ 6.2 ohm, tau ~ 0.25 ms.
        Self {
            turns: 500.0,
            radius_m: 0.01,
            length_m: 0.1,
            current_a: 2.0,
            frequency_hz: 1000.0,
            series_resistance_ohms: 4.0,
            core: CoreKind::Air,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Solenoid Coil Workbench right-side panel. A no-op when the
/// `show_coil_workbench` toggle is off.
pub fn draw_coil_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_coil_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_coil_workbench",
        "Solenoid Coil",
        |app, ui| {
            ui.label(
                egui::RichText::new("native long-solenoid inductor (field / L / X) · valenx-coil")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.coil;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Winding").strong());
                    ui.horizontal(|ui| {
                        ui.label("turns N");
                        ui.add(egui::DragValue::new(&mut s.turns).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("radius (m)");
                        ui.add(egui::DragValue::new(&mut s.radius_m).speed(0.001));
                    });
                    ui.horizontal(|ui| {
                        ui.label("length (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive").strong());
                    ui.horizontal(|ui| {
                        ui.label("current I (A)");
                        ui.add(egui::DragValue::new(&mut s.current_a).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_hz).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("series R (ohm)");
                        ui.add(egui::DragValue::new(&mut s.series_resistance_ohms).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Core").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.core, CoreKind::Air, "air");
                        ui.radio_value(&mut s.core, CoreKind::Iron, "iron");
                        ui.radio_value(&mut s.core, CoreKind::Ferrite, "ferrite");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_coil(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D coil").strong())
                        .on_hover_text(
                            "Build a representative solenoid winding (a cylinder shell with a thin inner core rod) as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Coil").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_coil_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.coil` borrow is released
    // here): build the coil's 3-D solid and load it.
    if app.coil.show_3d_request {
        app.coil.show_3d_request = false;
        load_coil_3d(app);
    }
}

/// Validate the form, evaluate the solenoid and format the readout.
fn run_coil(s: &mut CoilWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Solenoid`] from the form. The cross-sectional area
/// is `pi * r^2`. Extracted so it is unit-testable and shared with the 3-D
/// gate.
fn solenoid_of(s: &CoilWorkbenchState) -> Result<Solenoid, String> {
    let area = PI * s.radius_m * s.radius_m;
    Solenoid::new(s.turns, area, s.length_m, s.core.relative_permeability())
        .map_err(|e| e.to_string())
}

/// Evaluate the solenoid and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &CoilWorkbenchState) -> Result<String, String> {
    let coil = solenoid_of(s)?;
    let inductance = coil.inductance_henries();
    let energy = energy_joules(inductance, s.current_a);
    let reactance = reactance_ohms(inductance, s.frequency_hz).map_err(|e| e.to_string())?;
    // Series-RL time constant tau = L / R (the energising/decay e-folding
    // time): the current in a switched series-RL circuit reaches ~63.2 % of
    // its final value after one tau.
    let time_constant = coil
        .time_constant_seconds(s.series_resistance_ohms)
        .map_err(|e| e.to_string())?;
    // Axial field of an ideal long solenoid: B = mu0 * mu_r * (N / l) * I.
    let mu = VACUUM_PERMEABILITY * s.core.relative_permeability();
    let b_tesla = mu * (s.turns / s.length_m) * s.current_a;
    let ampere_turns = s.turns * s.current_a;
    let turns_per_metre = s.turns / s.length_m;

    Ok(format!(
        "turns N         : {:.0}\n\
         radius / length : {:.4} m / {:.4} m\n\
         core (mu_r)     : {} ({:.0})\n\
         current I       : {:.3} A\n\
         frequency f     : {:.1} Hz\n\
         series R        : {:.3} ohm\n\n\
         turn density n  : {:.1} turns/m\n\
         ampere-turns N·I: {:.2} A·turns\n\
         axial field B   : {:.6} T\n\
         inductance L    : {:.6} H\n\
         stored energy E : {:.6} J\n\
         reactance X_L   : {:.4} ohm\n\
         time const tau  : {:.6} s",
        s.turns,
        s.radius_m,
        s.length_m,
        s.core.label(),
        s.core.relative_permeability(),
        s.current_a,
        s.frequency_hz,
        s.series_resistance_ohms,
        turns_per_metre,
        ampere_turns,
        b_tesla,
        inductance,
        energy,
        reactance,
        time_constant,
    ))
}

/// Append a closed cylinder (axis along `z`, centre `c`, radius `r`,
/// half-height `hz`, `seg` radial segments) to the buffers as a double-sided
/// side wall plus two end caps.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hz: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Two rings of `seg` nodes (bottom then top), plus two cap centres.
    for ring in 0..2 {
        let z = if ring == 0 { c.z - hz } else { c.z + hz };
        for k in 0..seg {
            let a = TAU * (k as f64) / (seg as f64);
            nodes.push(Vector3::new(c.x + r * a.cos(), c.y + r * a.sin(), z));
        }
    }
    let bottom_centre = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, c.z - hz));
    let top_centre = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, c.z + hz));

    for k in 0..seg {
        let k1 = (k + 1) % seg;
        let b0 = base + k;
        let b1 = base + k1;
        let t0 = base + seg + k;
        let t1 = base + seg + k1;
        // Side wall, double-sided so it shades from either face.
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        tris.extend_from_slice(&[b0, t1, b1, b0, t0, t1]);
        // Caps.
        tris.extend_from_slice(&[bottom_centre, b1, b0]);
        tris.extend_from_slice(&[top_centre, t0, t1]);
    }
}

/// Build the solenoid as a triangle [`Mesh`] — the winding as a cylinder
/// shell (axis along `z`, length and radius from the form) with a thin inner
/// core rod. Representative geometry (true to the radius / length ratio; the
/// field / inductance numbers are the `valenx-coil` result). `None` for an
/// invalid configuration.
fn coil_solid_mesh(s: &CoilWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on the same validation the analysis uses.
    solenoid_of(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let half_len = 0.5 * s.length_m;
    let centre = Vector3::new(0.0, 0.0, half_len);

    // Winding cylinder.
    push_cylinder(&mut nodes, &mut tris, centre, s.radius_m, half_len, 28);
    // Thin inner core rod, slightly longer than the winding so it pokes out.
    push_cylinder(
        &mut nodes,
        &mut tris,
        centre,
        0.35 * s.radius_m,
        half_len * 1.15,
        20,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-coil");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D coil solid and load it into the central viewport.
fn load_coil_3d(app: &mut ValenxApp) {
    let Some(mesh) = coil_solid_mesh(&app.coil) else {
        app.coil.error = Some("coil parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<coil>/valenx-coil"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"coil"}`** product: the canonical
/// solenoid coil built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`CoilWorkbenchState::default`].
pub(crate) fn coil_product() -> crate::WorkspaceProduct {
    let s = CoilWorkbenchState::default();
    let mesh = coil_solid_mesh(&s).expect("canonical coil ⇒ solenoid solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<coil>/valenx-coil");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical coil ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Coil (solenoid)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = CoilWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_field_inductance_and_reactance() {
        let mut s = CoilWorkbenchState::default();
        run_coil(&mut s);
        assert!(
            s.error.is_none(),
            "default coil should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("axial field B"));
        assert!(s.result.contains("inductance L"));
        assert!(s.result.contains("reactance X_L"));
        assert!(s.result.contains("ampere-turns"));
        // 500 turns / 0.1 m at 2 A, air core: B = mu0 * 5000 * 2 ~ 12.6 mT.
        assert!(s.result.contains("0.012566"));
    }

    #[test]
    fn analyze_rejects_zero_length() {
        let mut s = CoilWorkbenchState {
            length_m: 0.0,
            ..Default::default()
        };
        run_coil(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn field_matches_hand_computed_air_core_solenoid() {
        // Ground truth: the axial field of an ideal long air-cored solenoid
        // is B = mu0 * N * I / l. For N = 1000, I = 3 A, l = 0.5 m:
        // B = (4 pi e-7) * 1000 * 3 / 0.5 = 4 pi e-7 * 6000 = 0.007539822 T.
        let s = CoilWorkbenchState {
            turns: 1000.0,
            radius_m: 0.02,
            length_m: 0.5,
            current_a: 3.0,
            core: CoreKind::Air,
            ..Default::default()
        };
        let coil = solenoid_of(&s).expect("valid air-core solenoid");
        // Inductance the crate computes, cross-checked: L = mu0 N^2 A / l.
        let area = PI * 0.02 * 0.02;
        let expected_l = VACUUM_PERMEABILITY * 1000.0 * 1000.0 * area / 0.5;
        assert!((coil.inductance_henries() - expected_l).abs() < 1.0e-12_f64 * expected_l);
        // Hand-computed axial field.
        let b: f64 = VACUUM_PERMEABILITY * 1000.0 * 3.0 / 0.5;
        let expected_b: f64 = 4.0 * PI * 1.0e-7 * 6000.0;
        assert!((b - expected_b).abs() < 1.0e-12);
    }

    #[test]
    fn time_constant_matches_hand_computed_and_is_in_readout() {
        // Ground truth: an air-cored solenoid with N = 100, A = 1e-4 m^2,
        // l = 0.1 m has L = mu0 * N^2 * A / l = 4 pi e-6 H. In series with
        // R = 2 ohm the RL time constant is tau = L / R = 2 pi e-6 s
        // = 6.283185e-6 s. Pick a radius giving A = pi r^2 = 1e-4 m^2, i.e.
        // r = sqrt(1e-4 / pi).
        let radius = (1.0e-4_f64 / PI).sqrt();
        let s = CoilWorkbenchState {
            turns: 100.0,
            radius_m: radius,
            length_m: 0.1,
            series_resistance_ohms: 2.0,
            core: CoreKind::Air,
            ..Default::default()
        };
        let coil = solenoid_of(&s).expect("valid air-core solenoid");
        let tau = coil
            .time_constant_seconds(s.series_resistance_ohms)
            .expect("positive R yields a time constant");
        let expected_tau: f64 = 2.0 * PI * 1.0e-6;
        assert!(
            (tau - expected_tau).abs() < 1.0e-15,
            "tau = {tau}, expected {expected_tau}"
        );
        // And the formatted readout surfaces it. The default coil
        // (N=500, r=0.01, l=0.1, air, R=4 ohm) has L ~ 0.987 mH and so
        // tau = L / R ~ 0.000247 s, formatted to 6 dp.
        let out = compute(&CoilWorkbenchState::default()).expect("default computes");
        assert!(out.contains("time const tau"), "readout: {out}");
        assert!(out.contains("0.000247"), "readout: {out}");
    }

    #[test]
    fn ferrite_core_multiplies_inductance_by_mu_r() {
        let air = CoilWorkbenchState {
            core: CoreKind::Air,
            ..Default::default()
        };
        let ferrite = CoilWorkbenchState {
            core: CoreKind::Ferrite,
            ..Default::default()
        };
        let l_air = solenoid_of(&air).unwrap().inductance_henries();
        let l_fe = solenoid_of(&ferrite).unwrap().inductance_henries();
        let ratio = l_fe / l_air;
        assert!((ratio - 2000.0).abs() < 1.0e-6, "ratio = {ratio}");
    }

    #[test]
    fn coil_mesh_for_default_is_nonempty_and_in_range() {
        let s = CoilWorkbenchState::default();
        let mesh = coil_solid_mesh(&s).expect("default coil yields a solid");
        assert!(mesh.nodes.len() > 8, "expected winding + core rings");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn coil_mesh_none_for_invalid() {
        let s = CoilWorkbenchState {
            length_m: 0.0,
            ..Default::default()
        };
        assert!(coil_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_coil_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_coil_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_coil_workbench = true;
        run_coil(&mut app.coil);
        draw_workbench(&mut app);
    }
}
