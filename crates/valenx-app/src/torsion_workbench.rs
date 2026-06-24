//! The right-side **Torsion Workbench** panel — native circular-shaft
//! elastic-torsion analysis over `valenx-torsion`.
//!
//! Mirrors the Heat Transfer / Flywheel workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_torsion_workbench`,
//! toggled from the View menu. The form picks a solid bar or a hollow tube,
//! sets the applied torque, length, shear modulus, spin speed and an
//! allowable shear stress, and "Analyze" evaluates the closed-form St.
//! Venant results — polar second moment `J`, the surface and at-radius shear
//! stress, the allowable torque, the angle of twist, the torsional rigidity
//! and the transmitted power. "Show 3-D" loads a representative shaft
//! cylinder (solid bar or bored tube) into the central viewport.

use std::f64::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_torsion::response::{
    allowable_torque, angle_of_twist, max_shear_stress, power, shear_stress_at, torsional_rigidity,
};
use valenx_torsion::shaft::Shaft;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which circular cross-section the form is describing.
#[derive(Copy, Clone, Debug, PartialEq)]
enum ShaftKind {
    /// A solid round bar of the outer diameter.
    Solid,
    /// A hollow round tube (a circular annulus) with a bore.
    Hollow,
}

/// Persistent form + result state for the Torsion Workbench.
pub struct TorsionWorkbenchState {
    /// Solid bar or hollow tube.
    kind: ShaftKind,
    /// Outer diameter `D` (m).
    outer_diameter_m: f64,
    /// Bore (inner) diameter `d` (m); used only for a hollow tube.
    inner_diameter_m: f64,
    /// Applied torque `T` (N·m).
    torque_nm: f64,
    /// Shaft length `L` over which the twist accumulates (m).
    length_m: f64,
    /// Shear modulus `G` of the material (Pa).
    shear_modulus_pa: f64,
    /// Angular speed `omega` for the power calculation (rad/s).
    angular_speed_rad_s: f64,
    /// Allowable shear stress for the torque-capacity check (Pa).
    allowable_shear_pa: f64,
    /// Radius at which to report the local shear stress (m); clamped into
    /// the cross-section before the query.
    query_radius_m: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D shaft solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for TorsionWorkbenchState {
    fn default() -> Self {
        // A 50 mm solid steel shaft (G = 79.3 GPa) carrying 1 kN·m at ~3000
        // rpm (omega ~ 314 rad/s), checked against a 60 MPa allowable: J =
        // pi d^4 / 32 ~ 6.136e-7 m^4, tau_max ~ 40.7 MPa, twist ~ 0.0021
        // rad/m, power ~ 314 kW.
        Self {
            kind: ShaftKind::Solid,
            outer_diameter_m: 0.05,
            inner_diameter_m: 0.03,
            torque_nm: 1000.0,
            length_m: 1.0,
            shear_modulus_pa: 79.3e9,
            angular_speed_rad_s: 314.159,
            allowable_shear_pa: 60.0e6,
            query_radius_m: 0.02,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Torsion Workbench right-side panel. A no-op when the
/// `show_torsion_workbench` toggle is off.
pub fn draw_torsion_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_torsion_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_torsion_workbench",
        "Torsion",
        |app, ui| {
            ui.label(
                egui::RichText::new("native circular-shaft elastic torsion · valenx-torsion")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.torsion;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Section").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.kind, ShaftKind::Solid, "solid bar");
                        ui.radio_value(&mut s.kind, ShaftKind::Hollow, "hollow tube");
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("outer diameter D (m)");
                        ui.add(egui::DragValue::new(&mut s.outer_diameter_m).speed(0.002))
                            .labelled_by(l.id);
                    });
                    if s.kind == ShaftKind::Hollow {
                        ui.horizontal(|ui| {
                            let l = ui.label("bore diameter d (m)");
                            ui.add(egui::DragValue::new(&mut s.inner_diameter_m).speed(0.002))
                                .labelled_by(l.id);
                        });
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loads").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("torque T (N·m)");
                        ui.add(egui::DragValue::new(&mut s.torque_nm).speed(10.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("length L (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(0.05))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("spin speed ω (rad/s)");
                        ui.add(egui::DragValue::new(&mut s.angular_speed_rad_s).speed(1.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material / limits").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("shear modulus G (Pa)");
                        ui.add(egui::DragValue::new(&mut s.shear_modulus_pa).speed(1.0e8))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("allowable τ (Pa)");
                        ui.add(egui::DragValue::new(&mut s.allowable_shear_pa).speed(1.0e6))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("query radius r (m)");
                        ui.add(egui::DragValue::new(&mut s.query_radius_m).speed(0.002))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_torsion(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative shaft cylinder (a solid bar, or a tube with its bore) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Torsion response").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_torsion_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.torsion` borrow is
    // released here): build the shaft's 3-D solid and load it.
    if app.torsion.show_3d_request {
        app.torsion.show_3d_request = false;
        load_shaft_3d(app);
    }
}

/// Validate the form, evaluate the shaft and format the readout.
fn run_torsion(s: &mut TorsionWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the [`Shaft`] for the current form, mapping any domain error to a
/// display string. Extracted so it is shared by the readout and the 3-D
/// gate.
fn build_shaft(s: &TorsionWorkbenchState) -> Result<Shaft, String> {
    match s.kind {
        ShaftKind::Solid => Shaft::solid(s.outer_diameter_m).map_err(|e| e.to_string()),
        ShaftKind::Hollow => {
            Shaft::hollow(s.outer_diameter_m, s.inner_diameter_m).map_err(|e| e.to_string())
        }
    }
}

/// Evaluate the shaft and format the full readout. Extracted so it is
/// unit-testable.
fn compute(s: &TorsionWorkbenchState) -> Result<String, String> {
    let shaft = build_shaft(s)?;
    let j = shaft.polar_moment();
    let zp = shaft.polar_section_modulus();
    let tau_max = max_shear_stress(&shaft, s.torque_nm).map_err(|e| e.to_string())?;

    // Clamp the query radius into the material before asking for the local
    // shear stress, so a stray entry reports the nearest valid stress
    // rather than an out-of-range error.
    let r_query = s
        .query_radius_m
        .clamp(shaft.inner_radius(), shaft.outer_radius());
    let tau_at = shear_stress_at(&shaft, s.torque_nm, r_query).map_err(|e| e.to_string())?;

    let t_allow = allowable_torque(&shaft, s.allowable_shear_pa).map_err(|e| e.to_string())?;
    let theta = angle_of_twist(&shaft, s.torque_nm, s.length_m, s.shear_modulus_pa)
        .map_err(|e| e.to_string())?;
    let theta_deg = theta * 180.0 / PI;
    let gj = torsional_rigidity(&shaft, s.shear_modulus_pa).map_err(|e| e.to_string())?;
    let p = power(s.torque_nm, s.angular_speed_rad_s).map_err(|e| e.to_string())?;

    let section = match s.kind {
        ShaftKind::Solid => format!("solid bar  D = {:.4} m", s.outer_diameter_m),
        ShaftKind::Hollow => format!(
            "hollow tube  D = {:.4} m  d = {:.4} m",
            s.outer_diameter_m, s.inner_diameter_m
        ),
    };

    Ok(format!(
        "{section}\n\
         torque T        : {:.1} N·m\n\
         length L        : {:.3} m\n\n\
         polar J         : {:.4e} m⁴\n\
         section mod Zp  : {:.4e} m³\n\
         tau_max (surface): {:.4e} Pa\n\
         tau at r={:.4} m : {:.4e} Pa\n\
         allowable torque: {:.1} N·m\n\
         angle of twist  : {:.6e} rad  ({:.4} deg)\n\
         rigidity G·J    : {:.4e} N·m²\n\
         power P=T·ω     : {:.4e} W",
        s.torque_nm, s.length_m, j, zp, tau_max, r_query, tau_at, t_allow, theta, theta_deg, gj, p,
    ))
}

/// Append a cylinder swept along the `z` axis: outer radius `r_out`, inner
/// (bore) radius `r_in` (`0` for a solid bar), spanning `z = 0 .. len`, to
/// the buffers as a closed triangulated solid (both end caps, the outer
/// wall, and the inner bore wall when `r_in > 0`). `seg` is the number of
/// angular segments.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    r_out: f64,
    r_in: f64,
    len: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Ring of nodes at each radius on the top (+z = len) and bottom (z = 0)
    // face. Order per angular step i: outer-top, outer-bottom, inner-top,
    // inner-bottom.
    for i in 0..seg {
        let a = std::f64::consts::TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(r_out * cos, r_out * sin, len));
        nodes.push(Vector3::new(r_out * cos, r_out * sin, 0.0));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, len));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, 0.0));
    }
    let mut quad = |a: usize, b: usize, c: usize, d: usize| {
        tris.extend_from_slice(&[base + a, base + b, base + c, base + a, base + c, base + d]);
    };
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (ot, ob, it, ib) = (4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
        let (not_, nob, nit, nib) = (4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3);
        // Top end-cap annulus, bottom end-cap annulus, outer wall, bore wall.
        quad(it, ot, not_, nit);
        quad(ob, ib, nib, nob);
        quad(ot, ob, nob, not_);
        quad(ib, it, nit, nib);
    }
}

/// Build the shaft as a triangle [`Mesh`] — a cylinder (solid bar, or a
/// tube with its bore) swept along its axis. Representative geometry: the
/// diameters are to scale, the length is a fixed visual aspect (the twist
/// and stress numbers are the `valenx-torsion` result). `None` for an
/// invalid configuration.
fn shaft_solid_mesh(s: &TorsionWorkbenchState) -> Option<Mesh> {
    let shaft = build_shaft(s).ok()?;
    let r_out = shaft.outer_radius();
    let r_in = shaft.inner_radius();
    // Representative axial length: a stubby aspect keyed to the diameter so
    // the section reads clearly in the viewport.
    let len = (8.0 * r_out).clamp(0.05, 2.0);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    push_cylinder(&mut nodes, &mut tris, r_out, r_in, len, 64);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-torsion");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D shaft solid and load it into the central viewport.
fn load_shaft_3d(app: &mut ValenxApp) {
    let Some(mesh) = shaft_solid_mesh(&app.torsion) else {
        app.torsion.error =
            Some("shaft parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<shaft>/valenx-torsion"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical torsion workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn torsion_product() -> crate::WorkspaceProduct {
    let s = TorsionWorkbenchState::default();
    let mesh = shaft_solid_mesh(&s).expect("canonical torsion ⇒ shaft solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<torsion>/valenx-shaft");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical torsion ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Shaft torsion (shear stress/twist)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
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

    #[test]
    fn default_state_is_idle() {
        let s = TorsionWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_every_quantity() {
        let mut s = TorsionWorkbenchState::default();
        run_torsion(&mut s);
        assert!(
            s.error.is_none(),
            "default shaft should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("polar J"));
        assert!(s.result.contains("tau_max (surface)"));
        assert!(s.result.contains("allowable torque"));
        assert!(s.result.contains("angle of twist"));
        assert!(s.result.contains("rigidity G·J"));
        assert!(s.result.contains("power P=T·ω"));
    }

    #[test]
    fn analyze_hollow_section_reports_tube_and_succeeds() {
        let mut s = TorsionWorkbenchState {
            kind: ShaftKind::Hollow,
            ..Default::default()
        };
        run_torsion(&mut s);
        assert!(
            s.error.is_none(),
            "hollow shaft should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("hollow tube"));
    }

    #[test]
    fn analyze_rejects_inverted_annulus() {
        // Bore not strictly smaller than the outer diameter is rejected by
        // the `valenx-torsion` constructor.
        let mut s = TorsionWorkbenchState {
            kind: ShaftKind::Hollow,
            outer_diameter_m: 0.03,
            inner_diameter_m: 0.05,
            ..Default::default()
        };
        run_torsion(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn surface_shear_stress_matches_t_r_over_j_solid() {
        // GROUND TRUTH: for a solid bar, tau_max = T * (D/2) / J with
        // J = pi * D^4 / 32, hand-computed independently of the crate.
        let d = 0.05_f64;
        let t = 1000.0_f64;
        let j = PI * d.powi(4) / 32.0;
        let expected = t * (d / 2.0) / j;

        let shaft = Shaft::solid(d).unwrap();
        let got = max_shear_stress(&shaft, t).unwrap();
        assert!(
            (got - expected).abs() < 1e-6 * expected,
            "tau_max got {got}, expected {expected}"
        );
    }

    #[test]
    fn query_radius_is_clamped_into_the_section() {
        // A query radius beyond the surface is clamped to the outer radius,
        // so the local stress equals the surface stress (and never errors).
        let s = TorsionWorkbenchState {
            query_radius_m: 10.0,
            ..Default::default()
        };
        let out = compute(&s).expect("clamped query should not error");
        assert!(out.contains("tau at r=0.0250 m"));
    }

    #[test]
    fn shaft_mesh_for_default_is_nonempty_and_in_range() {
        let s = TorsionWorkbenchState::default();
        let mesh = shaft_solid_mesh(&s).expect("default shaft yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a swept cylinder");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn shaft_mesh_none_for_invalid() {
        let s = TorsionWorkbenchState {
            outer_diameter_m: 0.0,
            ..Default::default()
        };
        assert!(shaft_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_torsion_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_torsion_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_torsion_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_torsion_workbench = true;
        run_torsion(&mut app.torsion);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every DragValue is a SpinButton that must be `labelled_by` its caption
        // (egui clears a DragValue's own Name), so an AI / screen reader can find
        // the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_torsion_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 7,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in [
            "outer diameter D (m)",
            "torque T (N·m)",
            "query radius r (m)",
        ] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
