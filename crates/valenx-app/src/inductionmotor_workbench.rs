//! The right-side **Induction Motor Workbench** panel — native 3-phase
//! induction-motor slip / power analysis over `valenx-inductionmotor`.
//!
//! Mirrors the Heat Pump / Battery Pack workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_inductionmotor_workbench`,
//! toggled from the View menu. The form sets the supply frequency, the
//! pole count, the rotor speed and the air-gap power; "Analyze" reports
//! the synchronous speed, slip, rotor frequency and the air-gap power
//! split into rotor copper loss and developed mechanical power, and
//! "Show 3-D motor" loads a representative TEFC motor solid into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;

use valenx_inductionmotor::InductionMotor;
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Blue-grey painted cast frame (a typical TEFC motor colour).
const FRAME: [f32; 3] = [0.22, 0.34, 0.52];
/// Grey machined end bells.
const ENDBELL: [f32; 3] = [0.52, 0.54, 0.57];
/// Dark cooling fins.
const FIN: [f32; 3] = [0.16, 0.18, 0.22];
/// Bright steel shaft + inner rotor.
const STEEL: [f32; 3] = [0.68, 0.70, 0.74];
/// Black plastic cooling-fan cowl.
const FAN: [f32; 3] = [0.10, 0.10, 0.11];
/// Black terminal box.
const TERMINAL: [f32; 3] = [0.10, 0.10, 0.12];
/// Dark mounting feet.
const FEET: [f32; 3] = [0.20, 0.21, 0.24];

/// Persistent form + result state for the Induction Motor Workbench.
pub struct InductionMotorWorkbenchState {
    /// Supply (line) frequency `f` (Hz).
    freq_hz: f64,
    /// Number of stator poles `p` (even, >= 2).
    poles: u32,
    /// Rotor mechanical speed `N` (rpm).
    rotor_rpm: f64,
    /// Air-gap power `P_airgap` crossing into the rotor (W).
    air_gap_power_w: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D motor solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for InductionMotorWorkbenchState {
    fn default() -> Self {
        // A 4-pole, 60 Hz induction motor at full-load slip (1750 rpm vs
        // 1800 sync = 2.78%), 10 kW across the air gap.
        Self {
            freq_hz: 60.0,
            poles: 4,
            rotor_rpm: 1750.0,
            air_gap_power_w: 10000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Induction Motor Workbench right-side panel. A no-op when the
/// `show_inductionmotor_workbench` toggle is off.
pub fn draw_inductionmotor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_inductionmotor_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_inductionmotor_workbench",
        "Induction Motor",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native 3-phase induction-motor slip / power · valenx-inductionmotor",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.inductionmotor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Supply + machine").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name otherwise, leaving it anonymous to a
                    // screen reader / AI driver).
                    ui.horizontal(|ui| {
                        let l = ui.label("frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.freq_hz).speed(0.5))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("poles (even)");
                        ui.add(egui::DragValue::new(&mut s.poles).speed(2.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("rotor speed N (rpm)");
                        ui.add(egui::DragValue::new(&mut s.rotor_rpm).speed(5.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("air-gap power (W)");
                        ui.add(egui::DragValue::new(&mut s.air_gap_power_w).speed(100.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_motor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D motor").strong())
                        .on_hover_text(
                            "Build a representative TEFC induction motor (frame, shaft, terminal box, fan cowl, feet) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_inductionmotor_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.inductionmotor` borrow
    // is released here): build the motor's 3-D solid and load it.
    if app.inductionmotor.show_3d_request {
        app.inductionmotor.show_3d_request = false;
        load_motor_3d(app);
    }
}

/// Validate the form, evaluate the machine and format the readout.
fn run_motor(s: &mut InductionMotorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a validated [`InductionMotor`] from the form, mapping the domain
/// error to a display string.
fn build_motor(s: &InductionMotorWorkbenchState) -> Result<InductionMotor, String> {
    InductionMotor::new(s.freq_hz, s.poles, s.rotor_rpm).map_err(|e| e.to_string())
}

/// Evaluate the machine and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &InductionMotorWorkbenchState) -> Result<String, String> {
    let m = build_motor(s)?;
    let cu_loss = m
        .rotor_copper_loss_w(s.air_gap_power_w)
        .map_err(|e| e.to_string())?;
    let mech = m
        .developed_mechanical_power_w(s.air_gap_power_w)
        .map_err(|e| e.to_string())?;
    Ok(format!(
        "supply          : {:.1} Hz, {} poles\n\
         sync speed Ns   : {:.1} rpm\n\
         rotor speed N   : {:.1} rpm\n\
         slip s          : {:.2} %  (slip speed {:.1} rpm)\n\
         rotor frequency : {:.2} Hz\n\n\
         air-gap power   : {:.2} kW\n\
         rotor Cu loss   : {:.2} kW\n\
         mech power Pmech: {:.2} kW",
        m.supply_frequency_hz(),
        m.poles(),
        m.sync_speed_rpm(),
        m.rotor_speed_rpm(),
        m.slip_percent(),
        m.slip_speed_rpm(),
        m.rotor_frequency_hz(),
        s.air_gap_power_w / 1000.0,
        cu_loss / 1000.0,
        mech / 1000.0,
    ))
}

/// Presentation spin rate of the rotor (inner rotor + shaft + cooling fan),
/// rad/s (~1.3 rev/s) — a readable inspect speed, not the real ~3000-rpm blur.
const ROTOR_RAD_PER_S: f32 = 8.0;

/// Build the induction motor as a triangle [`Mesh`] **with per-vertex colours**
/// plus the [`crate::RigidPart`] for the rotating rotor assembly. A
/// representative TEFC motor on the +x axle:
///
/// - blue-grey painted cast **frame** + a row of dark **cooling fins** along
///   the top;
/// - two grey machined **end bells**;
/// - the rotating **rotor assembly** — a visible steel **inner rotor** cylinder
///   inside the frame (so the spinning part reads, the fix for the old
///   monochrome/floating look), the protruding steel **output shaft** and the
///   black **fan cowl** — built consecutively as one contiguous node range that
///   spins about the motor axis;
/// - a black **terminal box** on top and dark mounting **feet**.
///
/// The frame / fins / end bells / terminal box / feet stay put while the rotor
/// spins. `None` for an invalid machine. Returns `(mesh, colors, parts)` with
/// `colors.len() == 3 × triangle_count`.
fn motor_solid_mesh_parts(
    s: &InductionMotorWorkbenchState,
) -> Option<(Mesh, Vec<[f32; 3]>, Vec<crate::RigidPart>)> {
    build_motor(s).ok()?;

    let axis = [1.0, 0.0, 0.0];
    let az = 0.5_f64; // axle height (z)
    let mut b = MeshBuilder::new();

    // Main frame (centre x = 0, length 1.4).
    b.cylinder([0.0, 0.0, az], axis, 0.42, 1.4, 28, FRAME);
    // Dark cooling fins: thin axial ribs along the top of the frame.
    for k in 0..6 {
        let fy = -0.18 + k as f64 * 0.072;
        b.cuboid([0.0, fy, az + 0.46], [1.2, 0.02, 0.12], FIN);
    }
    // End bells (slightly larger short cylinders at each end).
    b.cylinder([-0.73, 0.0, az], axis, 0.45, 0.1, 28, ENDBELL);
    b.cylinder([0.73, 0.0, az], axis, 0.45, 0.1, 28, ENDBELL);

    // Rotating rotor assembly: inner rotor (inside the frame) + output shaft +
    // fan cowl, built consecutively so they form one contiguous node range
    // spinning about the motor axis. Record its half-open span.
    let rotor_start = b.node_count();
    // Visible inner rotor cylinder inside the frame (the spinning part that was
    // previously implied, not shown).
    b.cylinder([0.0, 0.0, az], axis, 0.30, 1.3, 24, STEEL);
    // Output shaft (thin cylinder protruding from the +x end).
    b.cylinder([1.0, 0.0, az], axis, 0.08, 0.45, 16, STEEL);
    // Fan cowl at the −x end.
    b.cylinder([-0.86, 0.0, az], axis, 0.3, 0.17, 20, FAN);
    let rotor_end = b.node_count();

    // Terminal box on top (static).
    b.cuboid([0.0, 0.0, az + 0.46], [0.44, 0.36, 0.24], TERMINAL);
    // Mounting feet / base (static).
    b.cuboid([0.0, 0.0, 0.04], [1.4, 0.8, 0.08], FEET);

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-inductionmotor".to_string();

    // The rotor spins about the motor axis (+x) through the centreline
    // (y = 0, z = az).
    let parts = vec![crate::RigidPart {
        node_range: rotor_start..rotor_end,
        axis: [1.0, 0.0, 0.0],
        pivot: [0.0, 0.0, az as f32],
        rad_per_s: ROTOR_RAD_PER_S,
    }];
    Some((mesh, colors, parts))
}

/// Build the motor as a triangle [`Mesh`] (without the colour / rotor part
/// metadata) for the central viewport. See [`motor_solid_mesh_parts`].
fn motor_solid_mesh(s: &InductionMotorWorkbenchState) -> Option<Mesh> {
    motor_solid_mesh_parts(s).map(|(mesh, _colors, _parts)| mesh)
}

/// Build the 3-D motor solid and load it into the central viewport.
fn load_motor_3d(app: &mut ValenxApp) {
    let Some(mesh) = motor_solid_mesh(&app.inductionmotor) else {
        app.inductionmotor.error =
            Some("machine parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<motor>/valenx-inductionmotor"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical induction-motor workbench as a 3-D solid
/// plus its `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn inductionmotor_product() -> crate::WorkspaceProduct {
    let s = InductionMotorWorkbenchState::default();
    let (mesh, colors, parts) =
        motor_solid_mesh_parts(&s).expect("canonical induction motor ⇒ motor solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<inductionmotor>/valenx-motor");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical induction motor ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Induction motor (slip/torque)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: the rotor (output shaft + cooling fan) spins about the motor
        // axis while the frame/feet/terminal box stay put. Paused at t = 0.
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

    #[test]
    fn default_state_is_idle() {
        let s = InductionMotorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_slip_and_power() {
        let mut s = InductionMotorWorkbenchState::default();
        run_motor(&mut s);
        assert!(
            s.error.is_none(),
            "default motor should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("sync speed"));
        assert!(s.result.contains("slip"));
        assert!(s.result.contains("Pmech"));
        // 120 * 60 / 4 = 1800 rpm synchronous speed.
        assert!(s.result.contains("1800"));
    }

    #[test]
    fn analyze_rejects_odd_poles() {
        let mut s = InductionMotorWorkbenchState {
            poles: 3,
            ..Default::default()
        };
        run_motor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn air_gap_power_splits_into_copper_loss_and_mechanical() {
        // Ground truth: P_airgap = P_cu + P_mech with P_cu = s*P_ag and
        // P_mech = (1-s)*P_ag, and Ns = 120 f / p = 1800 rpm.
        let m = InductionMotor::new(60.0, 4, 1750.0).unwrap();
        assert!((m.sync_speed_rpm() - 1800.0).abs() < 1e-9);
        let pag = 10000.0;
        let cu = m.rotor_copper_loss_w(pag).unwrap();
        let mech = m.developed_mechanical_power_w(pag).unwrap();
        assert!(
            (cu + mech - pag).abs() < 1e-6,
            "cu {cu} + mech {mech} != {pag}"
        );
    }

    #[test]
    fn motor_mesh_for_default_is_nonempty_and_in_range() {
        let s = InductionMotorWorkbenchState::default();
        let mesh = motor_solid_mesh(&s).expect("default motor yields a solid");
        assert!(mesh.nodes.len() > 8, "expected frame + shaft + box + feet");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn motor_mesh_none_for_invalid() {
        let s = InductionMotorWorkbenchState {
            poles: 3,
            ..Default::default()
        };
        assert!(motor_solid_mesh(&s).is_none());
    }

    #[test]
    fn motor_carries_colours_and_a_visible_rotor() {
        // Per-vertex colours align to the renderer's coloured path (3/triangle).
        // The distinct part colours are present (frame / end bell / steel rotor /
        // fan / terminal / feet / fins) — the fix for the old monochrome look —
        // and the steel rotor assembly (inner rotor + shaft + fan) is a non-empty
        // mid-mesh range so the spinning part reads.
        let s = InductionMotorWorkbenchState::default();
        let (mesh, colors, parts) = motor_solid_mesh_parts(&s).expect("default motor builds");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        for col in [FRAME, ENDBELL, STEEL, FAN, TERMINAL, FEET, FIN] {
            assert!(colors.contains(&col), "missing part colour {col:?}");
        }
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
        // The rotor (inner rotor + shaft + fan) is a non-empty interior range.
        assert_eq!(parts.len(), 1);
        let p = &parts[0];
        assert!(p.node_range.start > 0 && p.node_range.end < mesh.nodes.len());
        assert_eq!(p.axis, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn inductionmotor_product_spins_the_rotor_only() {
        // The product carries a RigidParts animation: the rotor (shaft + fan, a
        // non-empty node range strictly inside the mesh) spins about +x at a
        // non-zero rate; the frame/feet/terminal box are left static.
        let product = inductionmotor_product();
        let loaded = product.mesh.as_ref().expect("motor product has a mesh");
        let node_count = loaded.mesh.nodes.len();
        let anim = product.animation.expect("motor product is animated");
        assert!(!anim.playing, "starts paused");
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => {
                assert_eq!(parts.len(), 1, "one rotating part: the rotor");
                let p = &parts[0];
                assert!(
                    p.node_range.start < p.node_range.end,
                    "non-empty rotor range"
                );
                assert!(
                    p.node_range.end <= node_count,
                    "rotor range within the mesh"
                );
                assert!(
                    p.node_range.start > 0 && p.node_range.end < node_count,
                    "frame precedes and feet/box follow the rotor (housing static)"
                );
                assert_eq!(p.axis, [1.0, 0.0, 0.0], "spins about the motor axis");
                assert!(p.rad_per_s.abs() > 0.0, "non-zero spin rate");
            }
            crate::ProductMotion::Turntable { .. } => {
                panic!("motor must use per-part rigid motion")
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
            draw_inductionmotor_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_inductionmotor_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_inductionmotor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_inductionmotor_workbench = true;
        run_motor(&mut app.inductionmotor);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Each numeric DragValue is a SpinButton; each must be `labelled_by`
        // its caption (egui clears a DragValue's own Name), so an AI / screen
        // reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_inductionmotor_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the inductionmotor numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every inductionmotor DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["frequency f (Hz)", "poles (even)", "rotor speed N (rpm)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
