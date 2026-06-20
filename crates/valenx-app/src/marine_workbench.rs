//! The right-side **Marine / Hull Workbench** panel — native box-form hull
//! hydrostatics + initial stability over `valenx-marine`.
//!
//! Mirrors the Springs workbench: a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_marine_workbench`, toggled from the View menu.
//! The form drives a [`valenx_marine::Hull`]; "Analyze" reports the
//! displacement, the centre of buoyancy `KB`, the transverse metacentric
//! radius `BM` and the metacentric height `GM` with a STABLE / UNSTABLE
//! verdict, and "Show 3-D hull" loads a raked-bow hull solid into the central
//! viewport (shaded, orbitable).

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_marine::{Hull, FRESHWATER_DENSITY, SEAWATER_DENSITY};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Marine / Hull Workbench.
pub struct MarineWorkbenchState {
    /// Waterline length `L` (m).
    length_m: f64,
    /// Beam / breadth `B` (m).
    beam_m: f64,
    /// Draft `T` (m).
    draft_m: f64,
    /// Block coefficient `Cb` in `(0, 1]`.
    block_coefficient: f64,
    /// Centre of gravity above the keel `KG` (m).
    kg_m: f64,
    /// Water density `rho` (kg/m^3).
    water_density: f64,
    /// Formatted hydrostatics readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D hull solid and load it into the
    /// central viewport (serviced after the panel draws).
    show_3d_request: bool,
}

impl Default for MarineWorkbenchState {
    fn default() -> Self {
        // A stable medium hull in seawater (GM > 0).
        Self {
            length_m: 120.0,
            beam_m: 20.0,
            draft_m: 6.0,
            block_coefficient: 0.70,
            kg_m: 8.0,
            water_density: SEAWATER_DENSITY,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Marine / Hull Workbench right-side panel. A no-op when the
/// `show_marine_workbench` toggle is off.
pub fn draw_marine_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_marine_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_marine_workbench",
        "Marine / Hull",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native box-form hull hydrostatics + stability · valenx-marine",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.marine;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Hull geometry (m)").strong());
                    ui.horizontal(|ui| {
                        ui.label("length L");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("beam B");
                        ui.add(egui::DragValue::new(&mut s.beam_m).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("draft T");
                        ui.add(egui::DragValue::new(&mut s.draft_m).speed(0.2));
                    });
                    ui.horizontal(|ui| {
                        ui.label("block coeff Cb");
                        ui.add(egui::DragValue::new(&mut s.block_coefficient).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loading + water").strong());
                    ui.horizontal(|ui| {
                        ui.label("KG (centre of gravity)");
                        ui.add(egui::DragValue::new(&mut s.kg_m).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("water ρ (kg/m³)");
                        ui.add(egui::DragValue::new(&mut s.water_density).speed(1.0));
                        if ui.small_button("sea").clicked() {
                            s.water_density = SEAWATER_DENSITY;
                        }
                        if ui.small_button("fresh").clicked() {
                            s.water_density = FRESHWATER_DENSITY;
                        }
                    });

                    // Live hint: block coefficient sweet spot by hull type.
                    ui.label(
                        egui::RichText::new(
                            "Cb ≈ 0.5 fine (yacht) · 0.7 cargo · 0.85+ tanker/barge",
                        )
                        .weak()
                        .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_marine(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D hull").strong())
                        .on_hover_text(
                            "Build the hull as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Hydrostatics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_marine_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.marine` borrow taken in
    // the closure above is released here): build the hull's 3-D solid and
    // load it into the central viewport, like the Springs workbench.
    if app.marine.show_3d_request {
        app.marine.show_3d_request = false;
        load_hull_3d(app);
    }
}

/// Build a validated [`Hull`] from the form, mapping the domain error to a
/// display string.
fn build_hull(s: &MarineWorkbenchState) -> Result<Hull, String> {
    Hull::new(
        s.length_m,
        s.beam_m,
        s.draft_m,
        s.block_coefficient,
        s.kg_m,
        s.water_density,
    )
    .map_err(|e| e.to_string())
}

/// Validate the form, compute the hydrostatics and format the readout.
/// Extracted from the draw closure so it is unit-testable.
fn run_marine(s: &mut MarineWorkbenchState) {
    s.error = None;
    match build_hull(s) {
        Ok(hull) => {
            let h = hull.hydrostatics();
            s.result = format!(
                "length L      : {:.2} m\n\
                 beam B        : {:.2} m\n\
                 draft T       : {:.2} m\n\
                 block coeff Cb: {:.3}\n\
                 KG            : {:.2} m\n\
                 water rho     : {:.0} kg/m³\n\n\
                 displaced vol : {:.1} m³\n\
                 displacement  : {:.1} t\n\
                 KB (buoyancy) : {:.3} m\n\
                 BM (metacentre): {:.3} m\n\
                 GM (metacentric height): {:.3} m\n\
                 stability     : {}",
                s.length_m,
                s.beam_m,
                s.draft_m,
                s.block_coefficient,
                s.kg_m,
                s.water_density,
                h.displaced_volume_m3,
                h.displacement_tonnes,
                h.kb_m,
                h.bm_m,
                h.gm_m,
                if h.stable {
                    "STABLE (GM > 0)"
                } else {
                    "UNSTABLE (GM <= 0)"
                },
            );
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build a representative ship hull as a triangle [`Mesh`] — a raked, pointed
/// bow forward (+x) tapering from the full midship / stern box section to a
/// stem at the waterline, with the keel forefoot set slightly aft of the
/// raked stem head. Length runs along x, beam along y, draft from the keel
/// `z = 0` up to the waterline `z = T`. Faces are emitted double-sided so the
/// shaded pass lights the hull from any orbit angle. The solid is a
/// representative hull *form*; the reported hydrostatics still use the
/// box-form `Cb` model from `valenx-marine`. `None` for an invalid hull.
fn hull_solid_mesh(s: &MarineWorkbenchState) -> Option<Mesh> {
    let hull = build_hull(s).ok()?;
    let (hl, hb, t) = (hull.length_m / 2.0, hull.beam_m / 2.0, hull.draft_m);
    // The full box section runs from the stern aft to `mx`; forward of that
    // the hull tapers in plan to the centreline bow. The keel forefoot
    // (`bkx`) sits aft of the raked stem head at the waterline (`hl`).
    let mx = hl * 0.15;
    let bkx = hl * 0.82;
    let nodes = vec![
        Vector3::new(-hl, -hb, 0.0), // 0 stern keel port
        Vector3::new(-hl, hb, 0.0),  // 1 stern keel stbd
        Vector3::new(-hl, hb, t),    // 2 stern deck stbd
        Vector3::new(-hl, -hb, t),   // 3 stern deck port
        Vector3::new(mx, -hb, 0.0),  // 4 mid keel port
        Vector3::new(mx, hb, 0.0),   // 5 mid keel stbd
        Vector3::new(mx, hb, t),     // 6 mid deck stbd
        Vector3::new(mx, -hb, t),    // 7 mid deck port
        Vector3::new(bkx, 0.0, 0.0), // 8 bow keel forefoot
        Vector3::new(hl, 0.0, t),    // 9 raked stem head (waterline)
    ];
    let mut tris: Vec<u32> = Vec::new();
    push_quad_ds(&mut tris, 0, 1, 2, 3); // transom (stern)
    push_quad_ds(&mut tris, 0, 1, 5, 4); // bottom, stern -> mid
    push_quad_ds(&mut tris, 3, 7, 6, 2); // deck, stern -> mid
    push_quad_ds(&mut tris, 0, 4, 7, 3); // port side, stern -> mid
    push_quad_ds(&mut tris, 1, 5, 6, 2); // stbd side, stern -> mid
    push_quad_ds(&mut tris, 4, 7, 9, 8); // port bow panel
    push_quad_ds(&mut tris, 5, 6, 9, 8); // stbd bow panel
    push_tri_ds(&mut tris, 4, 5, 8); // bottom forefoot wedge
    push_tri_ds(&mut tris, 7, 6, 9); // foredeck wedge to the stem
    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris;
    let mut mesh = Mesh::new("valenx-marine-hull");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Append a double-sided quad `a-b-c-d` (both windings) to `tris`.
fn push_quad_ds(tris: &mut Vec<u32>, a: usize, b: usize, c: usize, d: usize) {
    let (a, b, c, d) = (a as u32, b as u32, c as u32, d as u32);
    tris.extend_from_slice(&[a, b, c, a, c, d, a, c, b, a, d, c]);
}

/// Append a double-sided triangle `a-b-c` (both windings) to `tris`.
fn push_tri_ds(tris: &mut Vec<u32>, a: usize, b: usize, c: usize) {
    let (a, b, c) = (a as u32, b as u32, c as u32);
    tris.extend_from_slice(&[a, b, c, a, c, b]);
}

/// Build the 3-D hull solid and load it into the central viewport
/// (replacing any current STL / mesh) so it can be orbited — mirrors the
/// Springs workbench's `load_spring_3d`.
fn load_hull_3d(app: &mut ValenxApp) {
    let Some(mesh) = hull_solid_mesh(&app.marine) else {
        app.marine.error = Some("hull parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<hull>/valenx-marine"),
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
        let s = MarineWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_hull_is_stable() {
        let mut s = MarineWorkbenchState::default();
        run_marine(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        assert!(s.result.contains("displacement"));
        assert!(s.result.contains("GM"));
        assert!(s.result.contains("STABLE"));
    }

    #[test]
    fn analyze_rejects_bad_block_coefficient() {
        let mut s = MarineWorkbenchState {
            block_coefficient: 1.5,
            ..Default::default()
        };
        run_marine(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn a_high_centre_of_gravity_is_reported_unstable() {
        // KG well above the metacentre -> GM < 0.
        let mut s = MarineWorkbenchState {
            kg_m: 30.0,
            ..Default::default()
        };
        run_marine(&mut s);
        assert!(s.error.is_none());
        assert!(s.result.contains("UNSTABLE"));
    }

    #[test]
    fn hull_mesh_for_default_is_a_nonempty_hull() {
        let s = MarineWorkbenchState::default();
        let mesh = hull_solid_mesh(&s).expect("default hull yields a solid");
        assert_eq!(mesh.nodes.len(), 10);
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn hull_mesh_none_for_invalid_hull() {
        let s = MarineWorkbenchState {
            draft_m: 0.0,
            ..Default::default()
        };
        assert!(hull_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_marine_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_marine_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_marine_workbench = true;
        run_marine(&mut app.marine);
        draw_workbench(&mut app);
    }
}
