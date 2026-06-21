//! The right-side **Power Screw Workbench** panel — native square-thread
//! power / lead-screw analysis over `valenx-screwthread`.
//!
//! Mirrors the Torsion / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_screwthread_workbench`,
//! toggled from the View menu. The form sets the screw geometry (mean
//! diameter, pitch, thread starts), the screw/nut friction coefficient and an
//! axial load, and picks whether the load is being raised or lowered.
//! "Analyze" evaluates the closed-form square-thread power-screw results —
//! lead, lead (helix) angle, friction angle, the raising and lowering torque,
//! the thread-friction mechanical efficiency, and whether the screw is
//! self-locking. "Show 3-D" loads a representative screw-shank cylinder into
//! the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_screwthread::ScrewThread;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which direction the axial load is being driven, selecting the torque the
/// readout leads with.
#[derive(Copy, Clone, Debug, PartialEq)]
enum DriveMode {
    /// Drive the load against its direction (for example, lifting).
    Raise,
    /// Drive the load with its direction (for example, paying a weight down).
    Lower,
}

/// Persistent form + result state for the Power Screw Workbench.
pub struct ScrewThreadWorkbenchState {
    /// Mean (pitch-line) diameter `dm` (mm).
    mean_diameter_mm: f64,
    /// Single-thread pitch `p` (mm).
    pitch_mm: f64,
    /// Number of independent thread starts.
    starts: u32,
    /// Coefficient of friction `mu` between the screw and the nut.
    friction: f64,
    /// Axial load `F` resisted by the screw (N).
    axial_load_n: f64,
    /// Whether the load is being raised or lowered.
    mode: DriveMode,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D screw solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ScrewThreadWorkbenchState {
    fn default() -> Self {
        // Shigley single-start square-thread power screw lifting 6.5 kN:
        // dm = 25 mm, p = 5 mm, starts = 1, mu = 0.08. Self-locking, with a
        // raising torque of ~11732 N·mm and ~44% efficiency.
        Self {
            mean_diameter_mm: 25.0,
            pitch_mm: 5.0,
            starts: 1,
            friction: 0.08,
            axial_load_n: 6500.0,
            mode: DriveMode::Raise,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Power Screw Workbench right-side panel. A no-op when the
/// `show_screwthread_workbench` toggle is off.
pub fn draw_screwthread_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_screwthread_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_screwthread_workbench",
        "Power Screw",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native square-thread power / lead-screw torque · valenx-screwthread",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.screwthread;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Screw geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("mean diameter dm (mm)");
                        ui.add(egui::DragValue::new(&mut s.mean_diameter_mm).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("pitch p (mm)");
                        ui.add(egui::DragValue::new(&mut s.pitch_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("thread starts");
                        ui.add(egui::DragValue::new(&mut s.starts).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Friction & load").strong());
                    ui.horizontal(|ui| {
                        ui.label("friction μ");
                        ui.add(egui::DragValue::new(&mut s.friction).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("axial load F (N)");
                        ui.add(egui::DragValue::new(&mut s.axial_load_n).speed(10.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive direction").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, DriveMode::Raise, "raise");
                        ui.radio_value(&mut s.mode, DriveMode::Lower, "lower");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_screwthread(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative power-screw shank (a solid round bar of the mean diameter) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Power screw").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_screwthread_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.screwthread` borrow is
    // released here): build the screw's 3-D solid and load it.
    if app.screwthread.show_3d_request {
        app.screwthread.show_3d_request = false;
        load_screw_3d(app);
    }
}

/// Validate the form, evaluate the screw and format the readout.
fn run_screwthread(s: &mut ScrewThreadWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`ScrewThread`] from the form, mapping any domain
/// error to a display string. Extracted so it is shared and unit-testable.
fn build_screw(s: &ScrewThreadWorkbenchState) -> Result<ScrewThread, String> {
    ScrewThread::new(s.mean_diameter_mm, s.pitch_mm, s.starts, s.friction)
        .map_err(|e| e.to_string())
}

/// Evaluate the screw and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &ScrewThreadWorkbenchState) -> Result<String, String> {
    let screw = build_screw(s)?;
    let lead = screw.lead();
    let lead_angle_deg = screw.lead_angle().to_degrees();
    let friction_angle_deg = screw.friction_angle().to_degrees();
    let raise = screw
        .raise_torque(s.axial_load_n)
        .map_err(|e| e.to_string())?;
    let lower = screw.lower_torque(s.axial_load_n);
    let efficiency_pct = screw.efficiency().map_err(|e| e.to_string())? * 100.0;
    let self_locking = if screw.is_self_locking() { "yes" } else { "no" };
    let driving = match s.mode {
        DriveMode::Raise => raise,
        DriveMode::Lower => lower,
    };
    let mode_label = match s.mode {
        DriveMode::Raise => "raise",
        DriveMode::Lower => "lower",
    };

    Ok(format!(
        "mean diameter dm: {dm:.2} mm\n\
         pitch / starts  : {p:.2} mm / {starts}\n\
         friction μ      : {mu:.3}\n\
         axial load F    : {load:.1} N\n\n\
         lead l          : {lead:.3} mm\n\
         lead angle λ    : {lead_angle_deg:.4} °\n\
         friction angle φ: {friction_angle_deg:.4} °\n\
         self-locking    : {self_locking}\n\
         raise torque T_R: {raise:.2} N·mm\n\
         lower torque T_L: {lower:.2} N·mm\n\
         efficiency      : {efficiency_pct:.2} %\n\
         driving ({mode_label}) : {driving:.2} N·mm",
        dm = s.mean_diameter_mm,
        p = s.pitch_mm,
        starts = s.starts,
        mu = s.friction,
        load = s.axial_load_n,
    ))
}

/// Append a solid round bar (a closed z-cylinder) of radius `r` and length
/// `len`, spanning `z = 0..len`, to the buffers. `seg` angular segments.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    r: f64,
    len: f64,
    seg: usize,
) {
    let base = nodes.len();
    // Centre nodes for the two end-cap fans.
    let top_centre = base;
    let bot_centre = base + 1;
    nodes.push(Vector3::new(0.0, 0.0, len));
    nodes.push(Vector3::new(0.0, 0.0, 0.0));
    // Rim ring at each angular step: top then bottom.
    for i in 0..seg {
        let a = std::f64::consts::TAU * (i as f64) / (seg as f64);
        let (sin, cos) = a.sin_cos();
        nodes.push(Vector3::new(r * cos, r * sin, len));
        nodes.push(Vector3::new(r * cos, r * sin, 0.0));
    }
    let rim = base + 2;
    let mut tri = |a: usize, b: usize, c: usize| {
        tris.extend_from_slice(&[a, b, c]);
    };
    for i in 0..seg {
        let j = (i + 1) % seg;
        let (it, ib) = (rim + 2 * i, rim + 2 * i + 1);
        let (nit, nib) = (rim + 2 * j, rim + 2 * j + 1);
        // Top cap fan, bottom cap fan, then two wall triangles per segment.
        tri(top_centre, it, nit);
        tri(bot_centre, nib, ib);
        tri(it, ib, nib);
        tri(it, nib, nit);
    }
}

/// Build a representative power-screw shank as a triangle [`Mesh`] — a solid
/// round bar of the screw's mean diameter. Representative geometry (not to
/// scale; the torque / efficiency numbers are the `valenx-screwthread`
/// result). `None` for an invalid configuration.
fn screw_solid_mesh(s: &ScrewThreadWorkbenchState) -> Option<Mesh> {
    let screw = build_screw(s).ok()?;
    let r = screw.mean_diameter() / 2.0;
    // Representative axial length: several diameters of shank, clamped so the
    // bar reads clearly in the viewport whatever the diameter.
    let len = (6.0 * screw.mean_diameter()).clamp(20.0, 400.0);

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    push_cylinder(&mut nodes, &mut tris, r, len, 64);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-screwthread");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D screw solid and load it into the central viewport.
fn load_screw_3d(app: &mut ValenxApp) {
    let Some(mesh) = screw_solid_mesh(&app.screwthread) else {
        app.screwthread.error =
            Some("screw parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<screw>/valenx-screwthread"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"screwthread"}`** product: the
/// representative power-screw shank (a solid round bar of the mean diameter)
/// built from the canonical Shigley square-thread screw (dm = 25 mm, p = 5 mm,
/// single start, μ = 0.08, lifting 6.5 kN), paired with the lead-angle / torque
/// / efficiency / self-locking readout rows, at a fixed 3/4 camera. Registered
/// in [`crate::products_registry`]; the per-tool builder the registry
/// dispatches to. Pure — driven off [`ScrewThreadWorkbenchState::default`].
pub(crate) fn screwthread_product() -> crate::WorkspaceProduct {
    let s = ScrewThreadWorkbenchState::default();
    let mesh = screw_solid_mesh(&s).expect("canonical power screw ⇒ shank solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<screw>/valenx-screwthread");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical power screw ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Power screw (square thread)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = ScrewThreadWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
        assert_eq!(s.mode, DriveMode::Raise);
    }

    #[test]
    fn analyze_default_reports_torque_efficiency_and_self_locking() {
        let mut s = ScrewThreadWorkbenchState::default();
        run_screwthread(&mut s);
        assert!(
            s.error.is_none(),
            "default screw should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("raise torque T_R"));
        assert!(s.result.contains("lower torque T_L"));
        assert!(s.result.contains("efficiency"));
        // Shigley sample: T_R ~ 11732.29 N·mm, ~44.09% efficient, self-locks.
        assert!(s.result.contains("11732.29"));
        assert!(s.result.contains("44.09"));
        assert!(s.result.contains("self-locking    : yes"));
    }

    #[test]
    fn lower_mode_drives_with_the_lower_torque() {
        let mut s = ScrewThreadWorkbenchState {
            mode: DriveMode::Lower,
            ..Default::default()
        };
        run_screwthread(&mut s);
        assert!(s.error.is_none());
        // Lowering torque of the Shigley sample is ~1320.74 N·mm and is the
        // figure the "driving (lower)" line reports.
        assert!(s.result.contains("driving (lower) : 1320.74"));
    }

    #[test]
    fn analyze_rejects_zero_diameter() {
        let mut s = ScrewThreadWorkbenchState {
            mean_diameter_mm: 0.0,
            ..Default::default()
        };
        run_screwthread(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn lead_angle_matches_atan_ground_truth() {
        // Ground truth: lead = pitch * starts, and the lead (helix) angle is
        // atan(lead / (pi * dm)). Hand-computed for the default screw.
        let s = ScrewThreadWorkbenchState::default();
        let screw = build_screw(&s).expect("default screw is valid");
        let lead = s.pitch_mm * s.starts as f64;
        assert!((screw.lead() - lead).abs() < 1e-12);
        let expected_angle = (lead / (PI * s.mean_diameter_mm)).atan();
        assert!((screw.lead_angle() - expected_angle).abs() < 1e-12);
        // dm = 25, l = 5 -> tan(lambda) = 5 / (25*pi) = 0.0636620.
        let tan_lambda = (5.0_f64 / (25.0 * PI)).abs();
        assert!((screw.lead_angle().tan() - tan_lambda).abs() < 1e-9);
    }

    #[test]
    fn self_locking_exactly_when_friction_angle_exceeds_lead_angle() {
        // Ground truth: a screw self-locks iff its friction angle (atan mu)
        // is at least its lead angle. Default mu = 0.08 over a shallow helix
        // self-locks; a steep, low-friction multi-start screw back-drives.
        let locking = ScrewThreadWorkbenchState::default();
        let screw = build_screw(&locking).unwrap();
        assert!(screw.friction_angle() >= screw.lead_angle());
        assert!(screw.is_self_locking());

        let overhauling = ScrewThreadWorkbenchState {
            mean_diameter_mm: 20.0,
            pitch_mm: 12.0,
            starts: 4,
            friction: 0.02,
            ..Default::default()
        };
        let steep = build_screw(&overhauling).unwrap();
        assert!(steep.lead_angle() > steep.friction_angle());
        assert!(!steep.is_self_locking());
    }

    #[test]
    fn screw_mesh_for_default_is_nonempty_and_in_range() {
        let s = ScrewThreadWorkbenchState::default();
        let mesh = screw_solid_mesh(&s).expect("default screw yields a solid");
        assert!(mesh.nodes.len() > 8, "expected a faceted cylinder");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn screw_mesh_none_for_invalid() {
        let s = ScrewThreadWorkbenchState {
            pitch_mm: 0.0,
            ..Default::default()
        };
        assert!(screw_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_screwthread_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_screwthread_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_screwthread_workbench = true;
        run_screwthread(&mut app.screwthread);
        draw_workbench(&mut app);
    }
}
