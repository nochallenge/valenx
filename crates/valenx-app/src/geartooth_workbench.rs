//! The right-side **Gear Tooth Workbench** panel — native spur-gear tooth
//! bending-strength analysis over `valenx-geartooth`.
//!
//! Mirrors the Antenna / Heat-Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_geartooth_workbench`,
//! toggled from the View menu. The form sets the transmitted tangential
//! load, operating speed, face width, module and tooth count (plus the
//! AGMA quality `Qv`, fillet `Kf` and overload `Ko`); "Analyze" looks up
//! the Lewis form factor `Y(N)`, evaluates the Lewis root bending stress
//! `sigma = Wt/(F m Y)` and the AGMA refinement (whose dynamic factor `Kv`
//! is taken at the real pitch-line velocity `V = pi d n / 60000`), and
//! reports the pitch diameter, pitch-line velocity and both stresses, and
//! "Show 3-D gear" loads a representative gear-blank solid into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;

use valenx_gears::{GearKind, GearSpec};
use valenx_geartooth::agma::{agma_bending_stress, geometry_factor_j, AgmaFactors};
use valenx_geartooth::lewis::{lewis_bending_stress_for_teeth, pitch_line_velocity_m_per_s};
use valenx_geartooth::lewis_factor::lewis_form_factor;
use valenx_geartooth::spec::ToothLoad;
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Brushed-steel body colour for the gear blank.
const STEEL: [f32; 3] = [0.62, 0.64, 0.67];
/// A slightly darker cast colour for the centre hub / bore collar.
const HUB: [f32; 3] = [0.42, 0.44, 0.47];

/// Persistent form + result state for the Gear Tooth Workbench.
pub struct GeartoothWorkbenchState {
    /// Transmitted tangential load at the pitch line `Wt` (N).
    tangential_load_n: f64,
    /// Face width `F` (mm).
    face_width_mm: f64,
    /// Module `m` (mm). Pitch diameter = `module * teeth`.
    module_mm: f64,
    /// Number of teeth `N` on the gear (drives the Lewis form factor).
    teeth: u32,
    /// Pinion rotational speed `n` (rev/min). With the pitch diameter it
    /// fixes the pitch-line velocity `V = pi d n / 60000`, which feeds the
    /// AGMA dynamic factor `Kv`.
    speed_rpm: f64,
    /// AGMA transmission accuracy-level (quality) number `Qv` (6..=11).
    qv: f64,
    /// Fillet stress-concentration factor `Kf` (>= 1) for the AGMA `J`.
    fillet_kf: f64,
    /// AGMA overload factor `Ko` (>= 1).
    overload_ko: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D gear solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for GeartoothWorkbenchState {
    fn default() -> Self {
        // The crate's own doc example: a 20-tooth pinion, module 5 mm,
        // 50 mm face width, 3500 N tangential load. Shigley Table 14-2
        // gives Y = 0.322, so the Lewis stress is
        // 3500 / (50 * 5 * 0.322) ~= 43.48 MPa, and the pitch diameter is
        // module * teeth = 100 mm. Qv = 7 (commercial), Kf = 1.5, Ko = 1.25.
        // At 1000 rev/min the pitch-line velocity is
        // pi * 100 * 1000 / 60000 ~= 5.236 m/s.
        Self {
            tangential_load_n: 3500.0,
            face_width_mm: 50.0,
            module_mm: 5.0,
            teeth: 20,
            speed_rpm: 1000.0,
            qv: 7.0,
            fillet_kf: 1.5,
            overload_ko: 1.25,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Gear Tooth Workbench right-side panel. A no-op when the
/// `show_geartooth_workbench` toggle is off.
pub fn draw_geartooth_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_geartooth_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_geartooth_workbench",
        "Gear Tooth",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native spur-gear tooth bending strength (Lewis + AGMA) · valenx-geartooth",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.geartooth;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Gear geometry").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("module m (mm)");
                        ui.add(egui::DragValue::new(&mut s.module_mm).speed(0.1))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("teeth N");
                        ui.add(egui::DragValue::new(&mut s.teeth).speed(1.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("face width F (mm)");
                        ui.add(egui::DragValue::new(&mut s.face_width_mm).speed(0.5))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Operating load & speed").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("tangential Wt (N)");
                        ui.add(egui::DragValue::new(&mut s.tangential_load_n).speed(10.0))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("speed n (rev/min)");
                        ui.add(egui::DragValue::new(&mut s.speed_rpm).speed(10.0))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("AGMA factors").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("quality Qv (6..11)");
                        ui.add(egui::DragValue::new(&mut s.qv).speed(0.1))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("fillet Kf (≥1)");
                        ui.add(egui::DragValue::new(&mut s.fillet_kf).speed(0.05))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("overload Ko (≥1)");
                        ui.add(egui::DragValue::new(&mut s.overload_ko).speed(0.05))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_geartooth(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D gear").strong())
                        .on_hover_text(
                            "Build a representative spur-gear blank (disc, centre bore and a ring of teeth) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Tooth bending").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_geartooth_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.geartooth` borrow is
    // released here): build the gear's 3-D solid and load it.
    if app.geartooth.show_3d_request {
        app.geartooth.show_3d_request = false;
        load_gear_3d(app);
    }
}

/// Validate the form, evaluate the gear tooth and format the readout.
fn run_geartooth(s: &mut GeartoothWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the Lewis + AGMA bending stresses and format the full readout,
/// mapping any domain error to a display string. Extracted so it is
/// unit-testable.
fn compute(s: &GeartoothWorkbenchState) -> Result<String, String> {
    // Validate the load/geometry bundle (also gives the pitch diameter).
    let load = ToothLoad::new(s.tangential_load_n, s.face_width_mm, s.module_mm, s.teeth)
        .map_err(|e| e.to_string())?;
    let pitch_diameter_mm = load.pitch_diameter_mm();

    // Lewis: look up Y(N) and the bare bending stress.
    let lewis =
        lewis_bending_stress_for_teeth(s.tangential_load_n, s.face_width_mm, s.module_mm, s.teeth)
            .map_err(|e| e.to_string())?;

    // Pitch-line (tangential) velocity at the operating speed:
    // V = pi * d * n / 60000, in m/s. This is what the AGMA dynamic
    // factor curve-fit actually depends on.
    let velocity_m_per_s =
        pitch_line_velocity_m_per_s(pitch_diameter_mm, s.speed_rpm).map_err(|e| e.to_string())?;

    // AGMA refinement: J = Y / Kf, and the dynamic factor Kv from Qv at
    // the real pitch-line velocity V; here we apply the overload factor
    // Ko with the remaining factors at unity, replacing Y with J.
    let j = geometry_factor_j(lewis.form_factor_y, s.fillet_kf).map_err(|e| e.to_string())?;
    let kv = valenx_geartooth::agma::dynamic_factor_kv(s.qv, velocity_m_per_s)
        .map_err(|e| e.to_string())?;
    let factors = AgmaFactors::new(s.overload_ko, kv, 1.0, 1.0, 1.0).map_err(|e| e.to_string())?;
    let agma = agma_bending_stress(
        s.tangential_load_n,
        s.face_width_mm,
        s.module_mm,
        j,
        &factors,
    )
    .map_err(|e| e.to_string())?;

    Ok(format!(
        "module m        : {:.3} mm\n\
         teeth N         : {}\n\
         pitch diameter  : {:.3} mm\n\
         face width F    : {:.2} mm\n\
         tangential Wt   : {:.1} N\n\
         speed n         : {:.1} rev/min\n\
         pitch-line V    : {:.3} m/s\n\n\
         Lewis Y(N)      : {:.4}\n\
         Lewis sigma     : {:.3} MPa\n\n\
         fillet Kf       : {:.3}\n\
         geometry J=Y/Kf : {:.4}\n\
         dynamic Kv      : {:.4}\n\
         overload Ko     : {:.3}\n\
         AGMA sigma      : {:.3} MPa",
        s.module_mm,
        s.teeth,
        pitch_diameter_mm,
        s.face_width_mm,
        s.tangential_load_n,
        s.speed_rpm,
        velocity_m_per_s,
        lewis.form_factor_y,
        lewis.bending_stress_mpa,
        s.fillet_kf,
        j,
        kv,
        s.overload_ko,
        agma,
    ))
}

/// Build the spur gear as a triangle [`Mesh`] **with per-vertex colours** —
/// a real involute spur gear from [`valenx_gears::to_solid_spur`] (the same
/// quality bar as the Gears Workbench gear train), extruded along its axle
/// (+z) by the face width, plus a slightly darker centre hub/bore collar
/// ([`MeshBuilder::tube`]) so the bore reads. The gear geometry is the true
/// involute the bending-strength analysis is about — not the old ring of
/// rectangular bumps. Built in real millimetre units (the central viewport
/// auto-frames it); `None` for an invalid configuration.
///
/// Returns `(mesh, colors)` aligned to the renderer's coloured path
/// (`colors.len() == 3 × triangle_count`), so the product can set
/// [`crate::WorkspaceProduct::vertex_colors`].
fn gear_solid_mesh(s: &GeartoothWorkbenchState) -> Option<(Mesh, Vec<[f32; 3]>)> {
    // Gate on a valid load/geometry bundle (same check the readout uses) and a
    // tabulated Lewis form factor for this tooth count.
    let load = ToothLoad::new(s.tangential_load_n, s.face_width_mm, s.module_mm, s.teeth).ok()?;
    lewis_form_factor(s.teeth).ok()?;

    // The true involute spur gear (axle along +z, centred on z), tessellated.
    let spec = GearSpec {
        kind: GearKind::Spur,
        module_mm: s.module_mm,
        teeth: s.teeth,
        pressure_angle_deg: 20.0,
        helix_angle_deg: 0.0,
        face_width_mm: s.face_width_mm,
    };
    let solid = valenx_gears::to_solid_spur(&spec).ok()?;
    let gear_mesh = valenx_cad::solid_to_mesh(&solid, valenx_cad::DEFAULT_TESS_TOLERANCE).ok()?;

    // A darker hub/bore collar: a short annular tube on the same +z axle, bore
    // ~22 % of the pitch radius, standing slightly proud of each gear face so
    // it reads as a hub. The pitch radius is module·teeth/2 (mm).
    let pitch_r = load.pitch_diameter_mm() / 2.0;
    let bore_r = (pitch_r * 0.22).max(s.module_mm);
    let hub_r = bore_r + s.module_mm; // hub wall around the bore
    let hub_len = s.face_width_mm + 0.25 * s.face_width_mm;

    let mut b = MeshBuilder::new();
    b.append_tri_mesh(&gear_mesh, STEEL);
    b.tube(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        bore_r,
        hub_r,
        hub_len,
        28,
        HUB,
    );
    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-geartooth".to_string();
    Some((mesh, colors))
}

/// Build the 3-D gear solid and load it into the central viewport.
fn load_gear_3d(app: &mut ValenxApp) {
    let Some((mesh, _colors)) = gear_solid_mesh(&app.geartooth) else {
        app.geartooth.error =
            Some("gear parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<gear>/valenx-geartooth"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"geartooth"}`** product: the representative
/// spur-gear blank (disc + hub + tooth ring) built from the canonical
/// 20-tooth, module-5 mm, 3500 N example, paired with the Lewis + AGMA
/// tooth-bending readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to (so the geartooth 3-D product lives here, not in the shared reducer).
/// Pure (no app state) — driven entirely off [`GeartoothWorkbenchState::default`].
pub(crate) fn geartooth_product() -> crate::WorkspaceProduct {
    let s = GeartoothWorkbenchState::default();
    let (mesh, colors) =
        gear_solid_mesh(&s).expect("canonical gear-tooth spec ⇒ blank solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<gear>/valenx-geartooth");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical gear-tooth spec ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Spur gear tooth (Lewis + AGMA)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(colors),
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
        let s = GeartoothWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_stress_and_pitch_diameter() {
        let mut s = GeartoothWorkbenchState::default();
        run_geartooth(&mut s);
        assert!(
            s.error.is_none(),
            "default gear should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("pitch diameter"));
        assert!(s.result.contains("Lewis sigma"));
        assert!(s.result.contains("AGMA sigma"));
        // Default pitch diameter = module * teeth = 5 * 20 = 100 mm.
        assert!(s.result.contains("100.000 mm"));
        // Shigley Table 14-2: Y(20) = 0.322.
        assert!(s.result.contains("0.3220"));
    }

    #[test]
    fn analyze_rejects_zero_face_width() {
        let mut s = GeartoothWorkbenchState {
            face_width_mm: 0.0,
            ..Default::default()
        };
        run_geartooth(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn analyze_rejects_undercut_tooth_count() {
        // Below the minimum tabulated count (12), the Lewis factor is
        // out of domain, so the analyze must error.
        let mut s = GeartoothWorkbenchState {
            teeth: 8,
            ..Default::default()
        };
        run_geartooth(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_pitch_diameter_and_lewis_factor() {
        // Ground truth: pitch diameter = module * teeth, and the Lewis
        // form factor for a 20-tooth gear is exactly 0.322 (Shigley
        // Table 14-2), with sigma = Wt / (F m Y).
        let load = ToothLoad::new(3500.0, 50.0, 5.0, 20).unwrap();
        assert!((load.pitch_diameter_mm() - 100.0).abs() < 1e-9);

        let y = lewis_form_factor(20).unwrap();
        assert!((y - 0.322).abs() < 1e-9);

        let r = lewis_bending_stress_for_teeth(3500.0, 50.0, 5.0, 20).unwrap();
        let expected = 3500.0 / (50.0 * 5.0 * 0.322);
        assert!((r.bending_stress_mpa - expected).abs() < 1e-9);
    }

    #[test]
    fn analyze_default_reports_pitch_line_velocity() {
        use std::f64::consts::PI;

        let mut s = GeartoothWorkbenchState::default();
        run_geartooth(&mut s);
        assert!(
            s.error.is_none(),
            "default gear should analyze: {:?}",
            s.error
        );

        // The new readout lines: operating speed and pitch-line velocity.
        assert!(
            s.result.contains("speed n"),
            "missing speed line:\n{}",
            s.result
        );
        assert!(
            s.result.contains("pitch-line V"),
            "missing pitch-line velocity line:\n{}",
            s.result
        );
        // Default speed = 1000 rev/min (printed with one decimal place).
        assert!(
            s.result.contains("1000.0 rev/min"),
            "speed readout wrong:\n{}",
            s.result
        );

        // Ground truth: V = pi * d * n / 60000 with d = module * teeth =
        // 5 * 20 = 100 mm and n = 1000 rev/min, i.e.
        // pi * 100 * 1000 / 60000 = pi / 0.6 ~= 5.2359878 m/s, which the
        // {:.3} readout rounds to "5.236 m/s".
        let expected_v = PI * 100.0 * 1000.0 / 60_000.0;
        assert!(
            (expected_v - 5.235_987_755_982_99).abs() < 1e-9,
            "hand-computed V drifted: {expected_v}"
        );
        assert!(
            s.result.contains("5.236 m/s"),
            "pitch-line velocity readout wrong (expected ~{expected_v:.3} m/s):\n{}",
            s.result
        );
    }

    #[test]
    fn gear_mesh_for_default_is_nonempty_and_in_range() {
        let s = GeartoothWorkbenchState::default();
        let (mesh, _colors) = gear_solid_mesh(&s).expect("default gear yields a solid");
        assert!(mesh.nodes.len() > 8, "expected involute gear + hub");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn gear_mesh_carries_vertex_aligned_colours() {
        // The gear blank now ships per-vertex colours aligned to the renderer's
        // coloured path: exactly 3 per triangle (triangle-major). A drift here
        // would silently fall back to flat grey.
        let s = GeartoothWorkbenchState::default();
        let (mesh, colors) = gear_solid_mesh(&s).expect("default gear yields a solid");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
        // Two colours are present: the steel gear body and the darker hub.
        assert!(colors.contains(&STEEL), "gear body is steel");
        assert!(colors.contains(&HUB), "hub is a darker cast colour");
    }

    #[test]
    fn gear_mesh_none_for_invalid() {
        let s = GeartoothWorkbenchState {
            teeth: 0,
            ..Default::default()
        };
        assert!(gear_solid_mesh(&s).is_none());
    }

    #[test]
    fn geartooth_product_carries_colours() {
        let product = geartooth_product();
        let loaded = product.mesh.as_ref().expect("geartooth product has a mesh");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("geartooth product carries vertex_colors");
        assert_eq!(colors.len(), loaded.mesh.total_elements() * 3);
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
            draw_geartooth_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_geartooth_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_geartooth_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_geartooth_workbench = true;
        run_geartooth(&mut app.geartooth);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every DragValue is a SpinButton that must be `labelled_by` its caption
        // (egui clears a DragValue's own Name), so an AI / screen reader can find
        // the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_geartooth_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 8,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["module m (mm)", "tangential Wt (N)", "overload Ko (≥1)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
