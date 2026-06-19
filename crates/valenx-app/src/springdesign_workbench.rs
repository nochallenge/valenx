//! The right-side **Spring Design Workbench** panel — native round-wire
//! helical-compression-spring design over `valenx-spring-design`.
//!
//! Mirrors the Heat Transfer / BMR / Antenna workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_springdesign_workbench`,
//! toggled from the View menu. The form sets a spring's geometry (wire
//! diameter `d`, mean coil diameter `D`, active-coil count `N`), a material
//! whose preset fixes the shear modulus `G`, and an applied axial force;
//! "Analyze" builds a [`valenx_spring_design::HelicalSpring`] and reports
//! the spring index `C`, the Wahl curvature-correction factor, the spring
//! rate, the deflection under the force, and the Wahl-corrected (and
//! uncorrected) torsional shear stress, and "Show 3-D" loads a
//! representative helical coil solid into the central viewport to orbit.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_spring_design::HelicalSpring;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// A spring-wire material preset whose only role here is to fix the shear
/// modulus `G` (MPa) handed to [`HelicalSpring`]. Representative
/// engineering values for common round spring wire.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SpringMaterial {
    /// Music wire / spring steel, `G ~ 79_300 MPa`.
    SpringSteel,
    /// Austenitic stainless (302/304), `G ~ 69_000 MPa`.
    StainlessSteel,
    /// Phosphor bronze, `G ~ 41_400 MPa`.
    PhosphorBronze,
}

impl SpringMaterial {
    /// Shear modulus `G` (MPa) for this material preset.
    fn shear_modulus_mpa(self) -> f64 {
        match self {
            // Spring steel / music wire.
            SpringMaterial::SpringSteel => 79_300.0,
            // 302 / 304 austenitic stainless.
            SpringMaterial::StainlessSteel => 69_000.0,
            // Phosphor bronze (CuSn).
            SpringMaterial::PhosphorBronze => 41_400.0,
        }
    }

    /// Short human-readable label for the readout.
    fn label(self) -> &'static str {
        match self {
            SpringMaterial::SpringSteel => "spring steel",
            SpringMaterial::StainlessSteel => "stainless",
            SpringMaterial::PhosphorBronze => "phosphor bronze",
        }
    }
}

/// Persistent form + result state for the Spring Design Workbench.
pub struct SpringDesignWorkbenchState {
    /// Wire diameter `d` (mm).
    wire_diameter_mm: f64,
    /// Mean coil diameter `D`, centre-to-centre of the wire (mm).
    mean_coil_diameter_mm: f64,
    /// Number of *active* (deflecting) coils `N`.
    active_coils: f64,
    /// Material preset, which fixes the shear modulus `G`.
    material: SpringMaterial,
    /// Applied axial force `F` (N).
    force_n: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D coil solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for SpringDesignWorkbenchState {
    fn default() -> Self {
        // A 2 mm music-wire spring, 16 mm mean coil, 10 active coils under
        // 40 N: spring index C = 8, rate ~ 3.87 N/mm, deflection ~ 10.3 mm,
        // Wahl-corrected shear stress ~ 241 MPa.
        Self {
            wire_diameter_mm: 2.0,
            mean_coil_diameter_mm: 16.0,
            active_coils: 10.0,
            material: SpringMaterial::SpringSteel,
            force_n: 40.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Spring Design Workbench right-side panel. A no-op when the
/// `show_springdesign_workbench` toggle is off.
pub fn draw_springdesign_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_springdesign_workbench {
        return;
    }

    egui::SidePanel::right("valenx_springdesign_workbench")
        .resizable(true)
        .default_width(360.0)
        .width_range(300.0..=560.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Spring Design",
                "native helical compression spring design · valenx-spring-design",
            ) {
                app.show_springdesign_workbench = false;
            }

            let s = &mut app.springdesign;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        ui.label("wire d (mm)");
                        ui.add(egui::DragValue::new(&mut s.wire_diameter_mm).speed(0.05));
                    });
                    ui.horizontal(|ui| {
                        ui.label("mean coil D (mm)");
                        ui.add(egui::DragValue::new(&mut s.mean_coil_diameter_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("active coils N");
                        ui.add(egui::DragValue::new(&mut s.active_coils).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Material (sets G)").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(
                            &mut s.material,
                            SpringMaterial::SpringSteel,
                            "spring steel",
                        );
                        ui.radio_value(
                            &mut s.material,
                            SpringMaterial::StainlessSteel,
                            "stainless",
                        );
                    });
                    ui.radio_value(
                        &mut s.material,
                        SpringMaterial::PhosphorBronze,
                        "phosphor bronze",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Load").strong());
                    ui.horizontal(|ui| {
                        ui.label("applied force F (N)");
                        ui.add(egui::DragValue::new(&mut s.force_n).speed(1.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_springdesign(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative helical compression coil as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Spring performance").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        });

    // Serviced after the panel draws (the `&mut app.springdesign` borrow is
    // released here): build the coil's 3-D solid and load it.
    if app.springdesign.show_3d_request {
        app.springdesign.show_3d_request = false;
        load_spring_3d(app);
    }
}

/// Validate the form, evaluate the spring and format the readout.
fn run_springdesign(s: &mut SpringDesignWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`HelicalSpring`] for the current form, the value
/// both the readout and the 3-D gate need. Extracted so it is unit-testable
/// and shared.
fn build_spring(s: &SpringDesignWorkbenchState) -> Result<HelicalSpring, String> {
    HelicalSpring::new(
        s.wire_diameter_mm,
        s.mean_coil_diameter_mm,
        s.active_coils,
        s.material.shear_modulus_mpa(),
    )
    .map_err(|e| e.to_string())
}

/// Evaluate the spring and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &SpringDesignWorkbenchState) -> Result<String, String> {
    let spring = build_spring(s)?;
    let c = spring.spring_index();
    let kw = spring.wahl_factor();
    let k = spring.rate();
    let deflection = spring.deflection(s.force_n).map_err(|e| e.to_string())?;
    let tau = spring.shear_stress(s.force_n).map_err(|e| e.to_string())?;
    let tau0 = spring
        .shear_stress_uncorrected(s.force_n)
        .map_err(|e| e.to_string())?;

    // Advisory note when the index is outside the conventional 4..=12 band.
    let index_note = if spring.index_is_typical() {
        "typical (4..12)"
    } else {
        "outside 4..12 band"
    };

    Ok(format!(
        "wire d / coil D : {:.3} / {:.3} mm\n\
         active coils N  : {:.2}\n\
         material (G)    : {} ({:.0} MPa)\n\
         applied force F : {:.2} N\n\n\
         spring index C  : {:.3}  [{index_note}]\n\
         Wahl factor Kw  : {:.4}\n\
         spring rate k   : {:.4} N/mm\n\
         deflection delta: {:.3} mm\n\
         shear tau (Wahl): {:.2} MPa\n\
         shear tau0 (raw): {:.2} MPa",
        s.wire_diameter_mm,
        s.mean_coil_diameter_mm,
        s.active_coils,
        s.material.label(),
        s.material.shear_modulus_mpa(),
        s.force_n,
        c,
        kw,
        k,
        deflection,
        tau,
        tau0,
    ))
}

/// Append an axis-aligned box (centre `c`, half-extents `h`) to the
/// buffers, as 12 outward-facing triangles.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Build the helical compression spring as a triangle [`Mesh`] — the wire
/// swept along a helix of mean radius `D/2` over `N` turns, approximated by
/// a chain of small wire-thick boxes, climbing in `z`. Representative
/// geometry (the engineering numbers are the `valenx-spring-design`
/// result, not this mesh). `None` for an invalid configuration.
fn coil_solid_mesh(s: &SpringDesignWorkbenchState) -> Option<Mesh> {
    // Reuse the crate's validation so the mesh and the readout agree on
    // what counts as a constructible spring.
    build_spring(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Work in mm directly; the viewport auto-frames the result.
    let radius = s.mean_coil_diameter_mm / 2.0;
    let wire = s.wire_diameter_mm;
    // A solid coil needs the pitch >= the wire diameter so turns do not
    // overlap; use 1.5*d so the helix reads clearly when orbited.
    let pitch = (1.5 * wire).max(wire);
    // Sample the helix densely enough that the box chain looks smooth, and
    // clamp the total so a huge coil count cannot explode the buffer.
    let turns = s.active_coils;
    let segments = ((turns * 24.0).round() as usize).clamp(24, 4096);
    // Half-extents of each wire segment box: half the wire thickness in the
    // cross-section, and a little extra along the path so adjacent boxes
    // overlap into a continuous tube.
    let seg_half = wire * 0.6;
    let r_half = wire * 0.5;

    for i in 0..=segments {
        let t = i as f64 / segments as f64; // 0..1 over the whole coil
        let angle = t * turns * std::f64::consts::TAU;
        let x = radius * angle.cos();
        let y = radius * angle.sin();
        let z = t * turns * pitch;
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(x, y, z),
            Vector3::new(r_half, r_half, seg_half),
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-spring-design");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D coil solid and load it into the central viewport.
fn load_spring_3d(app: &mut ValenxApp) {
    let Some(mesh) = coil_solid_mesh(&app.springdesign) else {
        app.springdesign.error =
            Some("spring parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<coil>/valenx-spring-design"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = SpringDesignWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_index_rate_and_stress() {
        let mut s = SpringDesignWorkbenchState::default();
        run_springdesign(&mut s);
        assert!(
            s.error.is_none(),
            "default spring should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("spring index C"));
        assert!(s.result.contains("Wahl factor Kw"));
        assert!(s.result.contains("spring rate k"));
        // d=2, D=16 -> C = 8.000; k = 3.8721 N/mm; deflection = 10.330 mm;
        // Wahl-corrected shear = 241.21 MPa.
        assert!(s.result.contains("8.000"));
        assert!(s.result.contains("3.8721"));
        assert!(s.result.contains("10.330"));
        assert!(s.result.contains("241.21"));
    }

    #[test]
    fn analyze_rejects_wire_thicker_than_coil() {
        // D < d -> degenerate spring (index <= 1).
        let mut s = SpringDesignWorkbenchState {
            wire_diameter_mm: 20.0,
            mean_coil_diameter_mm: 16.0,
            ..Default::default()
        };
        run_springdesign(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn analyze_rejects_zero_force() {
        let mut s = SpringDesignWorkbenchState {
            force_n: 0.0,
            ..Default::default()
        };
        run_springdesign(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn spring_index_and_wahl_match_hand_computed_ground_truth() {
        // Ground truth: spring index C = D/d, and the Wahl curvature factor
        // Kw = (4C - 1)/(4C - 4) + 0.615/C. For d=2, D=16: C = 8 and
        // Kw = 31/28 + 0.615/8 = 1.1840178571...
        let s = SpringDesignWorkbenchState::default();
        let spring = build_spring(&s).expect("default spring is valid");
        assert!((spring.spring_index() - 8.0).abs() < 1e-9);
        let c = 16.0_f64 / 2.0_f64;
        let kw_hand = (4.0 * c - 1.0) / (4.0 * c - 4.0) + 0.615 / c;
        assert!((spring.wahl_factor() - kw_hand).abs() < 1e-12);
        assert!((spring.wahl_factor() - 1.1840178571428).abs() < 1e-9);
    }

    #[test]
    fn coil_mesh_for_default_is_nonempty_and_in_range() {
        let s = SpringDesignWorkbenchState::default();
        let mesh = coil_solid_mesh(&s).expect("default spring yields a coil solid");
        assert!(mesh.nodes.len() > 8, "expected a multi-segment coil");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn coil_mesh_none_for_invalid() {
        let s = SpringDesignWorkbenchState {
            wire_diameter_mm: 20.0,
            mean_coil_diameter_mm: 16.0,
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
            draw_springdesign_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_springdesign_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_springdesign_workbench = true;
        run_springdesign(&mut app.springdesign);
        draw_workbench(&mut app);
    }
}
