//! The right-side **Flywheel Workbench** panel — native closed-form
//! flywheel (rotational kinetic energy storage) sizing over
//! `valenx-flywheel`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_flywheel_workbench`,
//! toggled from the View menu. The form sets a rotor (shape, mass, radii)
//! of a chosen material density spinning between an upper and a lower
//! operating speed (rpm); "Analyze" reports the moment of inertia, the
//! stored kinetic energy at the top speed, the usable energy delivered as
//! the rotor slows to the bottom speed, the coefficient of fluctuation of
//! that speed band, and the rim speed and first-order hoop stress at the
//! rim; "Show 3-D" loads a representative spinning disc (solid or annular)
//! solid into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;
use std::f64::consts::TAU;

use valenx_flywheel::{coefficient_of_fluctuation, rpm_to_rad_s, Flywheel, Rotor};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// The rotor geometry chosen in the form.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum RotorShape {
    /// Solid uniform disk / cylinder, `I = 1/2 m r^2`.
    SolidDisk,
    /// Thin ring / rim with all mass at the radius, `I = m r^2`.
    ThinRing,
    /// Hollow annular disk between `r_in` and `r_out`,
    /// `I = 1/2 m (r_in^2 + r_out^2)`.
    AnnularDisk,
}

/// Persistent form + result state for the Flywheel Workbench.
pub struct FlywheelWorkbenchState {
    /// The rotor geometry.
    shape: RotorShape,
    /// Rotor mass `m` (kg).
    mass_kg: f64,
    /// Outer radius `r_out` (m). Also the rim radius for a ring / disk.
    r_out_m: f64,
    /// Inner radius `r_in` (m), used by the annular disk only.
    r_in_m: f64,
    /// Material density `rho` (kg/m^3), used by the rim-stress model.
    density_kg_m3: f64,
    /// Upper operating speed (rev/min).
    rpm_max: f64,
    /// Lower operating speed (rev/min).
    rpm_min: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D rotor solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for FlywheelWorkbenchState {
    fn default() -> Self {
        // A 10 kg solid steel disk, 0.3 m radius (rho = 7800 kg/m^3),
        // cycling 3000 -> 1500 rpm. I = 0.45 kg.m^2; at 3000 rpm
        // (~314 rad/s) it stores ~22 kJ and delivers ~17 kJ down to
        // 1500 rpm, with a wide coefficient of fluctuation of ~0.67.
        Self {
            shape: RotorShape::SolidDisk,
            mass_kg: 10.0,
            r_out_m: 0.3,
            r_in_m: 0.15,
            density_kg_m3: 7800.0,
            rpm_max: 3000.0,
            rpm_min: 1500.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Flywheel Workbench right-side panel. A no-op when the
/// `show_flywheel_workbench` toggle is off.
pub fn draw_flywheel_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_flywheel_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_flywheel_workbench",
        "Flywheel",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form flywheel energy + rim-stress sizing · valenx-flywheel",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.flywheel;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Rotor shape").strong());
                    ui.radio_value(&mut s.shape, RotorShape::SolidDisk, "solid disk / cylinder");
                    ui.radio_value(&mut s.shape, RotorShape::ThinRing, "thin ring / rim");
                    ui.radio_value(&mut s.shape, RotorShape::AnnularDisk, "annular (hollow) disk");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise).
                    ui.horizontal(|ui| {
                        let mk = ui.label("mass (kg)");
                        ui.add(egui::DragValue::new(&mut s.mass_kg).speed(0.1))
                            .labelled_by(mk.id)
                            .on_hover_text("Flywheel mass (kg)");
                    });
                    ui.horizontal(|ui| {
                        let ro = ui.label("outer radius (m)");
                        ui.add(egui::DragValue::new(&mut s.r_out_m).speed(0.01))
                            .labelled_by(ro.id)
                            .on_hover_text("Outer radius (m)");
                    });
                    if s.shape == RotorShape::AnnularDisk {
                        ui.horizontal(|ui| {
                            let ri = ui.label("inner radius (m)");
                            ui.add(egui::DragValue::new(&mut s.r_in_m).speed(0.01))
                                .labelled_by(ri.id)
                                .on_hover_text("Inner (bore) radius (m)");
                        });
                    }
                    ui.horizontal(|ui| {
                        let de = ui.label("density ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.density_kg_m3).speed(10.0))
                            .labelled_by(de.id)
                            .on_hover_text("Material density ρ (kg/m³)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating band").strong());
                    ui.horizontal(|ui| {
                        let rx = ui.label("upper speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rpm_max).speed(10.0))
                            .labelled_by(rx.id)
                            .on_hover_text("Upper operating speed (rpm)");
                    });
                    ui.horizontal(|ui| {
                        let rn = ui.label("lower speed (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rpm_min).speed(10.0))
                            .labelled_by(rn.id)
                            .on_hover_text("Lower operating speed (rpm)");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_flywheel(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build the rotor as a representative spinning disc solid (solid or annular) and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Energy + stress").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_flywheel_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.flywheel` borrow is
    // released here): build the rotor's 3-D solid and load it.
    if app.flywheel.show_3d_request {
        app.flywheel.show_3d_request = false;
        load_flywheel_3d(app);
    }
}

/// Validate the form, evaluate the flywheel and format the readout.
fn run_flywheel(s: &mut FlywheelWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Flywheel`] for the current form, the object both
/// the readout and the 3-D gate need. Extracted so it is unit-testable and
/// shared.
fn build_flywheel(s: &FlywheelWorkbenchState) -> Result<Flywheel, String> {
    let rotor = match s.shape {
        RotorShape::SolidDisk => Rotor::solid_disk(s.mass_kg, s.r_out_m),
        RotorShape::ThinRing => Rotor::thin_ring(s.mass_kg, s.r_out_m),
        RotorShape::AnnularDisk => Rotor::annular_disk(s.mass_kg, s.r_in_m, s.r_out_m),
    }
    .map_err(|e| e.to_string())?;
    Flywheel::new(rotor, s.density_kg_m3).map_err(|e| e.to_string())
}

/// Evaluate the flywheel and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &FlywheelWorkbenchState) -> Result<String, String> {
    let fw = build_flywheel(s)?;
    let omega_max = rpm_to_rad_s(s.rpm_max).map_err(|e| e.to_string())?;
    let omega_min = rpm_to_rad_s(s.rpm_min).map_err(|e| e.to_string())?;

    let inertia = fw.moment_of_inertia();
    let energy_top = fw.energy_at(omega_max).map_err(|e| e.to_string())?;
    // Specific energy E/m at the top speed — the energy-density figure of
    // merit (J/kg) that lets a rotor's storage be compared independent of
    // how heavy it is.
    let specific = fw.specific_energy(omega_max).map_err(|e| e.to_string())?;
    let usable = fw
        .usable_energy(omega_min, omega_max)
        .map_err(|e| e.to_string())?;
    let cs = coefficient_of_fluctuation(omega_min, omega_max).map_err(|e| e.to_string())?;
    let v_rim = fw.rim_speed(omega_max).map_err(|e| e.to_string())?;
    let sigma = fw.rim_stress(omega_max).map_err(|e| e.to_string())?;

    let shape = match s.shape {
        RotorShape::SolidDisk => "solid disk",
        RotorShape::ThinRing => "thin ring",
        RotorShape::AnnularDisk => "annular disk",
    };

    Ok(format!(
        "rotor shape     : {shape}\n\
         mass            : {:.3} kg\n\
         outer radius    : {:.4} m\n\
         density ρ       : {:.1} kg/m³\n\
         speed band      : {:.0} → {:.0} rpm\n\n\
         inertia I       : {inertia:.5} kg·m²\n\
         energy @ top    : {:.1} J\n\
         specific energy : {specific:.1} J/kg\n\
         usable energy   : {:.1} J\n\
         coeff. of fluct.: {cs:.4}\n\
         rim speed       : {v_rim:.2} m/s\n\
         rim hoop stress : {:.3} MPa",
        s.mass_kg,
        s.r_out_m,
        s.density_kg_m3,
        s.rpm_max,
        s.rpm_min,
        energy_top,
        usable,
        sigma / 1.0e6,
    ))
}

/// Append a swept disc — a flat ring of outer radius `r_out` and inner
/// radius `r_in` (0 for a solid disk), centred on the spin axis (the `z`
/// axis) and of axial half-thickness `hz` — to the buffers as a closed
/// triangulated solid (top, bottom, outer wall, and the inner bore wall
/// when `r_in > 0`). `seg` is the number of angular segments.
fn push_disc(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    r_out: f64,
    r_in: f64,
    hz: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Ring of nodes at each radius, on the top (+z) and bottom (-z) face.
    // Order per angular step i: outer-top, outer-bottom, inner-top,
    // inner-bottom.
    for i in 0..seg {
        let a = TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(r_out * cos, r_out * sin, hz));
        nodes.push(Vector3::new(r_out * cos, r_out * sin, -hz));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, hz));
        nodes.push(Vector3::new(r_in * cos, r_in * sin, -hz));
    }
    let mut quad = |a: usize, b: usize, c: usize, d: usize| {
        tris.extend_from_slice(&[base + a, base + b, base + c, base + a, base + c, base + d]);
    };
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (ot, ob, it, ib) = (4 * i, 4 * i + 1, 4 * i + 2, 4 * i + 3);
        let (not_, nob, nit, nib) = (4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3);
        // Top annulus, bottom annulus, outer wall, inner bore wall.
        quad(it, ot, not_, nit);
        quad(ob, ib, nib, nob);
        quad(ot, ob, nob, not_);
        quad(ib, it, nit, nib);
    }
}

/// Build the rotor as a triangle [`Mesh`] — a spinning disc (solid for the
/// disk / ring shapes, annular for the hollow disk) about the `z` axis.
/// Representative geometry (a fixed visual thickness; the energy and stress
/// numbers are the `valenx-flywheel` result). `None` for an invalid
/// configuration.
/// Presentation spin rate of the flywheel disc, rad/s (~1.3 rev/s) — a readable
/// inspect speed, not the real ~3000-rpm blur.
const ROTOR_RAD_PER_S: f32 = 8.0;

/// Build the rotor as a triangle [`Mesh`] together with the
/// [`crate::RigidPart`] for the spinning disc. The whole mesh *is* the rotor (a
/// single body), so the one part covers every node, spinning about the +z axis
/// through the origin. `None` for an invalid configuration.
fn rotor_solid_mesh_parts(s: &FlywheelWorkbenchState) -> Option<(Mesh, Vec<crate::RigidPart>)> {
    build_flywheel(s).ok()?;

    let r_out = s.r_out_m;
    let r_in = if s.shape == RotorShape::AnnularDisk {
        s.r_in_m
    } else {
        0.0
    };
    // Representative axial half-thickness: a slim disc keyed to the radius.
    let hz = (0.08 * r_out).clamp(0.005, 0.05);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    push_disc(&mut nodes, &mut tris, r_out, r_in, hz, 64);

    let node_count = nodes.len();
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-flywheel");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    // The whole disc is the rotor: one part covering every node, spinning about
    // the disc (+z) axis through the origin.
    let parts = vec![crate::RigidPart {
        node_range: 0..node_count,
        axis: [0.0, 0.0, 1.0],
        pivot: [0.0, 0.0, 0.0],
        rad_per_s: ROTOR_RAD_PER_S,
    }];
    Some((mesh, parts))
}

/// Build the rotor as a triangle [`Mesh`] (without the rotor part metadata) for
/// the central viewport. See [`rotor_solid_mesh_parts`].
fn rotor_solid_mesh(s: &FlywheelWorkbenchState) -> Option<Mesh> {
    rotor_solid_mesh_parts(s).map(|(mesh, _parts)| mesh)
}

/// Build the 3-D rotor solid and load it into the central viewport.
fn load_flywheel_3d(app: &mut ValenxApp) {
    let Some(mesh) = rotor_solid_mesh(&app.flywheel) else {
        app.flywheel.error =
            Some("rotor parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<rotor>/valenx-flywheel"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"flywheel"}`** product: the representative
/// spinning rotor disc built from the canonical 10 kg / 0.3 m solid steel disk
/// cycling 3000→1500 rpm, paired with the energy + rim-stress readout rows, at
/// a fixed 3/4 camera. Registered in [`crate::products_registry`]; the
/// per-tool builder the registry dispatches to. Pure — driven off
/// [`FlywheelWorkbenchState::default`].
pub(crate) fn flywheel_product() -> crate::WorkspaceProduct {
    let s = FlywheelWorkbenchState::default();
    let (mesh, parts) =
        rotor_solid_mesh_parts(&s).expect("canonical flywheel ⇒ rotor disc solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<rotor>/valenx-flywheel");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical flywheel ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Flywheel (energy + rim stress)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: the disc spins about its own (+z) axis. A single rigid part
        // covering the whole rotor (set explicitly so the default turntable does
        // not apply). Paused at t = 0.
        animation: Some(crate::ProductAnimation {
            playing: false,
            speed: 1.0,
            t: 0.0,
            motion: crate::ProductMotion::RigidParts(parts),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_flywheel::rad_s_to_rpm;

    #[test]
    fn default_state_is_idle() {
        let s = FlywheelWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_inertia_energy_and_stress() {
        let mut s = FlywheelWorkbenchState::default();
        run_flywheel(&mut s);
        assert!(
            s.error.is_none(),
            "default rotor should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("inertia I"));
        assert!(s.result.contains("usable energy"));
        assert!(s.result.contains("rim hoop stress"));
        // 10 kg, 0.3 m solid disk: I = 0.5 * 10 * 0.09 = 0.45 kg.m^2.
        assert!(s.result.contains("0.45000"));
    }

    #[test]
    fn analyze_rejects_zero_mass() {
        let mut s = FlywheelWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        run_flywheel(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_inner_ge_outer_annulus() {
        let mut s = FlywheelWorkbenchState {
            shape: RotorShape::AnnularDisk,
            r_in_m: 0.4,
            r_out_m: 0.3,
            ..Default::default()
        };
        run_flywheel(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn kinetic_energy_is_half_i_omega_squared_ground_truth() {
        // Ground truth, hand-computed: a 10 kg / 0.3 m solid disk has
        // I = 1/2 m r^2 = 0.45 kg.m^2; at 3000 rpm,
        // omega = 3000 * 2 pi / 60 = 100 pi rad/s, so
        // E = 1/2 I omega^2 = 0.5 * 0.45 * (100 pi)^2.
        let s = FlywheelWorkbenchState::default();
        let fw = build_flywheel(&s).unwrap();
        let inertia: f64 = fw.moment_of_inertia();
        assert!((inertia - 0.45).abs() < 1e-12);
        let omega: f64 = rpm_to_rad_s(3000.0).unwrap();
        let expected: f64 = 0.5 * 0.45 * omega * omega;
        let energy = fw.energy_at(omega).unwrap();
        assert!((energy - expected).abs() < 1e-6, "got {energy}");
        // Round-trip the speed conversion as an independent check.
        let back: f64 = rad_s_to_rpm(omega).unwrap();
        assert!((back - 3000.0).abs() < 1e-9);
    }

    #[test]
    fn specific_energy_is_energy_over_mass_ground_truth() {
        // GROUND TRUTH, hand-computed: the 10 kg / 0.3 m solid disk stores
        // E = 1/2 I omega^2 with I = 0.45 kg.m^2 at omega = 3000 rpm
        // (= 100 pi rad/s), so the energy density is
        // E/m = (0.5 * 0.45 * (100 pi)^2) / 10 = 2250 pi^2 / 10 J/kg.
        let s = FlywheelWorkbenchState::default();
        let fw = build_flywheel(&s).unwrap();
        let omega: f64 = rpm_to_rad_s(3000.0).unwrap();
        let expected: f64 = (0.5 * 0.45 * omega * omega) / 10.0;
        let se = fw.specific_energy(omega).unwrap();
        assert!((se - expected).abs() < 1e-6, "specific energy got {se}");
        // Pi-pinned closed form, independent of the crate's arithmetic path.
        let two_pi = 2.0_f64 * std::f64::consts::PI;
        let analytic: f64 = 2250.0 * (two_pi / 2.0).powi(2) / 10.0;
        assert!((se - analytic).abs() < 1e-6, "analytic se got {se}");

        // The readout surfaces it at one decimal: 2250 pi^2 / 10 = 2220.7.
        let mut shown = FlywheelWorkbenchState::default();
        run_flywheel(&mut shown);
        assert!(
            shown.result.contains("specific energy : 2220.7 J/kg"),
            "readout was: {}",
            shown.result
        );
    }

    #[test]
    fn rotor_mesh_for_default_is_nonempty_and_in_range() {
        let s = FlywheelWorkbenchState::default();
        let mesh = rotor_solid_mesh(&s).expect("default rotor yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a swept disc");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn rotor_mesh_none_for_invalid() {
        let s = FlywheelWorkbenchState {
            mass_kg: 0.0,
            ..Default::default()
        };
        assert!(rotor_solid_mesh(&s).is_none());
    }

    #[test]
    fn flywheel_product_spins_the_disc() {
        // The product carries a RigidParts animation: the whole disc (one part
        // covering every node) spins about the +z disc axis at a non-zero rate.
        let product = flywheel_product();
        let loaded = product.mesh.as_ref().expect("flywheel product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("flywheel product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 1, "one rotating part: the disc");
                let p = &parts[0];
                assert_eq!(p.node_range.start, 0, "the disc is the whole mesh");
                assert_eq!(p.node_range.end, node_count, "range covers every node");
                assert_eq!(p.axis, [0.0, 0.0, 1.0], "spins about the disc axis");
                assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("flywheel must set RigidParts explicitly so the default no-ops")
            }
        }
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
            draw_flywheel_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_flywheel_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_flywheel_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_flywheel_workbench = true;
        run_flywheel(&mut app.flywheel);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The geometry + operating-band DragValues are SpinButtons; each must
        // be `labelled_by` its caption (egui clears a DragValue's own Name), so
        // an AI / screen reader can find the control by the caption text. The
        // default shape is a solid disk, so the inner-radius control is hidden;
        // mass, outer radius, density, upper speed and lower speed remain.
        let mut app = ValenxApp::default();
        app.show_flywheel_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // mass, outer radius, density, upper speed, lower speed.
        assert!(
            spin_buttons.len() >= 5,
            "expected the flywheel numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every flywheel DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["mass (kg)", "outer radius (m)", "upper speed (rpm)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
