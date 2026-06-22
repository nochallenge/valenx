//! The right-side **Insulation Workbench** panel — native composite-wall
//! U-value and R-value analysis over `valenx-insulation`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_insulation_workbench`,
//! toggled from the View menu. The form describes a two-layer building wall
//! (an insulation layer plus a structural layer) sandwiched between an
//! interior and exterior ISO-6946 surface film; "Analyze" sums the
//! area-specific series resistances and reports the per-layer R-values, the
//! total R-value, the overall U-value and the steady-state heat-loss rate,
//! and "Show 3-D" loads a representative layered wall slab into the central
//! viewport.
//!
//! This wraps the `valenx-insulation` building-envelope model, whose scope
//! differs from the `valenx-heat-transfer` workbench: here the films are the
//! fixed ISO-6946 reference resistances (`R_si = 0.13`, `R_se = 0.04`
//! `m^2.K/W`) and the readout is framed around the area-specific R-value /
//! U-value an insulation specifier works in, rather than convective film
//! coefficients.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_insulation::{CompositeWall, Layer, SurfaceFilm};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Insulation Workbench.
pub struct InsulationWorkbenchState {
    /// Insulation-layer thickness `L1` (m).
    insulation_thickness_m: f64,
    /// Insulation-layer conductivity `k1` (W/m·K).
    insulation_k_w_per_mk: f64,
    /// Structural-layer thickness `L2` (m).
    structure_thickness_m: f64,
    /// Structural-layer conductivity `k2` (W/m·K).
    structure_k_w_per_mk: f64,
    /// Include the ISO-6946 interior surface film (`R_si = 0.13`).
    include_interior_film: bool,
    /// Include the ISO-6946 exterior surface film (`R_se = 0.04`).
    include_exterior_film: bool,
    /// Wall area `A` (m^2), used for the heat-loss rate.
    area_m2: f64,
    /// Inside-minus-outside temperature difference `dT` (K).
    delta_t_k: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D wall solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for InsulationWorkbenchState {
    fn default() -> Self {
        // A 120 mm k=0.035 insulation layer plus a 100 mm k=0.15 block,
        // between the ISO-6946 reference films, over 10 m^2 at dT = 20 K:
        // R ~ 4.27 m^2.K/W, U ~ 0.23 W/m^2K, ~47 W of loss.
        Self {
            insulation_thickness_m: 0.12,
            insulation_k_w_per_mk: 0.035,
            structure_thickness_m: 0.10,
            structure_k_w_per_mk: 0.15,
            include_interior_film: true,
            include_exterior_film: true,
            area_m2: 10.0,
            delta_t_k: 20.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Insulation Workbench right-side panel. A no-op when the
/// `show_insulation_workbench` toggle is off.
pub fn draw_insulation_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_insulation_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_insulation_workbench",
        "Insulation",
        |app, ui| {
            ui.label(
                egui::RichText::new("native composite-wall R-value / U-value · valenx-insulation")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.insulation;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Insulation layer").strong());
                    ui.horizontal(|ui| {
                        ui.label("thickness (m)");
                        ui.add(egui::DragValue::new(&mut s.insulation_thickness_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("conductivity k (W/m·K)");
                        ui.add(egui::DragValue::new(&mut s.insulation_k_w_per_mk).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Structural layer").strong());
                    ui.horizontal(|ui| {
                        ui.label("thickness (m)");
                        ui.add(egui::DragValue::new(&mut s.structure_thickness_m).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("conductivity k (W/m·K)");
                        ui.add(egui::DragValue::new(&mut s.structure_k_w_per_mk).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Surface films (ISO 6946)").strong());
                    ui.checkbox(&mut s.include_interior_film, "interior film (R_si 0.13)");
                    ui.checkbox(&mut s.include_exterior_film, "exterior film (R_se 0.04)");

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Heat loss").strong());
                    ui.horizontal(|ui| {
                        ui.label("area (m²)");
                        ui.add(egui::DragValue::new(&mut s.area_m2).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ΔT inside−outside (K)");
                        ui.add(egui::DragValue::new(&mut s.delta_t_k).speed(0.5));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_insulation(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative layered wall slab (insulation + structural layer, with the ISO-6946 surface films) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("U-value").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_insulation_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.insulation` borrow is
    // released here): build the wall's 3-D solid and load it.
    if app.insulation.show_3d_request {
        app.insulation.show_3d_request = false;
        load_wall_3d(app);
    }
}

/// Validate the form, evaluate the wall and format the readout.
fn run_insulation(s: &mut InsulationWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Assemble the [`CompositeWall`] for the current form — an insulation
/// layer, a structural layer, and the optional ISO-6946 interior/exterior
/// films. Extracted so the readout and the 3-D gate share one builder.
fn build_wall(s: &InsulationWorkbenchState) -> Result<CompositeWall, String> {
    let insulation =
        Layer::new(s.insulation_thickness_m, s.insulation_k_w_per_mk).map_err(|e| e.to_string())?;
    let structure =
        Layer::new(s.structure_thickness_m, s.structure_k_w_per_mk).map_err(|e| e.to_string())?;

    let mut builder = CompositeWall::builder();
    if s.include_interior_film {
        builder = builder.interior_film(SurfaceFilm::interior_default());
    }
    builder = builder.layer(insulation).layer(structure);
    if s.include_exterior_film {
        builder = builder.exterior_film(SurfaceFilm::exterior_default());
    }
    builder.build().map_err(|e| e.to_string())
}

/// Evaluate the wall and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &InsulationWorkbenchState) -> Result<String, String> {
    let wall = build_wall(s)?;

    let r_insulation = wall.layers()[0].resistance();
    let r_structure = wall.layers()[1].resistance();
    let r_total = wall.total_resistance();
    let u = wall.u_value();
    let q = wall
        .heat_loss(s.area_m2, s.delta_t_k)
        .map_err(|e| e.to_string())?;

    let r_si = wall.interior_film().map(|f| f.resistance()).unwrap_or(0.0);
    let r_se = wall.exterior_film().map(|f| f.resistance()).unwrap_or(0.0);

    Ok(format!(
        "insulation L / k: {:.3} m / {:.3} W/m·K\n\
         structure  L / k: {:.3} m / {:.3} W/m·K\n\
         area / ΔT       : {:.2} m² / {:.1} K\n\n\
         R interior film : {:.4} m²K/W\n\
         R insulation    : {:.4} m²K/W\n\
         R structure     : {:.4} m²K/W\n\
         R exterior film : {:.4} m²K/W\n\
         R total         : {:.4} m²K/W\n\
         U-value         : {:.3} W/m²K\n\
         heat loss Q     : {:.2} W",
        s.insulation_thickness_m,
        s.insulation_k_w_per_mk,
        s.structure_thickness_m,
        s.structure_k_w_per_mk,
        s.area_m2,
        s.delta_t_k,
        r_si,
        r_insulation,
        r_structure,
        r_se,
        r_total,
        u,
        q,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
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

/// Build the composite wall as a triangle [`Mesh`] — a layered slab (the
/// insulation layer plus the structural layer stacked through the wall's `x`
/// direction), with thin inside / outside surface films and a base.
/// Representative geometry (not to scale; the R-value / U-value numbers are
/// the `valenx-insulation` result). `None` for an invalid configuration.
fn wall_solid_mesh(s: &InsulationWorkbenchState) -> Option<Mesh> {
    build_wall(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Insulation layer (inner half of the slab, -x side).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.05, 0.0, 0.6),
        Vector3::new(0.05, 0.7, 0.5),
    );
    // Structural layer (outer half of the slab, +x side).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.05, 0.0, 0.6),
        Vector3::new(0.05, 0.7, 0.5),
    );
    // Interior surface film (-x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.115, 0.0, 0.6),
        Vector3::new(0.012, 0.66, 0.46),
    );
    // Exterior surface film (+x face).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.115, 0.0, 0.6),
        Vector3::new(0.012, 0.66, 0.46),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.25, 0.7, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-insulation");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D wall solid and load it into the central viewport.
fn load_wall_3d(app: &mut ValenxApp) {
    let Some(mesh) = wall_solid_mesh(&app.insulation) else {
        app.insulation.error =
            Some("wall parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<wall>/valenx-insulation"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: a DATA-ONLY text card of the insulation workbench's
/// `compute()` readout rows (see [`crate::products_registry`]). An R-value /
/// U-value result has no characteristic shape — the panel's layered slab is a
/// fixed schematic stack of boxes (not scaled to the real layer thicknesses),
/// not a real object — so the bridge product is right-sized to a card
/// (`mesh: None`) carrying just the readout (the confidence badge is appended
/// centrally). The panel's "Show 3-D" button still builds that representative
/// slab into the central viewport.
pub(crate) fn insulation_product() -> crate::WorkspaceProduct {
    let s = InsulationWorkbenchState::default();
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical insulation ⇒ readout computes"),
    );
    crate::WorkspaceProduct {
        title: "Insulation (R-value/heat loss)".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: Default::default(),
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
        let s = InsulationWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_rvalue_uvalue_and_heat_loss() {
        let mut s = InsulationWorkbenchState::default();
        run_insulation(&mut s);
        assert!(
            s.error.is_none(),
            "default wall should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("R total"));
        assert!(s.result.contains("U-value"));
        assert!(s.result.contains("heat loss Q"));
        // 120 mm of k=0.035 dominates: U ~ 0.23 W/m^2K.
        assert!(s.result.contains("0.23"));
    }

    #[test]
    fn analyze_rejects_zero_conductivity() {
        let mut s = InsulationWorkbenchState {
            insulation_k_w_per_mk: 0.0,
            ..Default::default()
        };
        run_insulation(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn u_value_is_reciprocal_of_total_and_q_is_u_a_dt() {
        // Ground truth, hand-computed: insulation R = 0.12/0.035, structure
        // R = 0.10/0.15, plus R_si = 0.13 and R_se = 0.04; U = 1/R_total and
        // Q = U * A * dT exactly.
        let s = InsulationWorkbenchState::default();
        let wall = build_wall(&s).unwrap();
        let r_total = 0.13 + 0.12 / 0.035 + 0.10 / 0.15 + 0.04;
        assert!((wall.total_resistance() - r_total).abs() < 1e-9);
        let u: f64 = 1.0 / r_total;
        assert!((wall.u_value() - u).abs() < 1e-9);
        assert!((wall.u_value() * wall.total_resistance() - 1.0).abs() < 1e-12);
        let q = wall.heat_loss(s.area_m2, s.delta_t_k).unwrap();
        assert!((q - u * s.area_m2 * s.delta_t_k).abs() < 1e-9);
    }

    #[test]
    fn dropping_films_lowers_total_resistance() {
        // Removing the surface films removes R_si + R_se = 0.17 m^2.K/W.
        let with_films = InsulationWorkbenchState::default();
        let no_films = InsulationWorkbenchState {
            include_interior_film: false,
            include_exterior_film: false,
            ..Default::default()
        };
        let r_with = build_wall(&with_films).unwrap().total_resistance();
        let r_without = build_wall(&no_films).unwrap().total_resistance();
        assert!(r_with > r_without);
        assert!((r_with - r_without - 0.17).abs() < 1e-9);
    }

    #[test]
    fn wall_mesh_for_default_is_nonempty_and_in_range() {
        let s = InsulationWorkbenchState::default();
        let mesh = wall_solid_mesh(&s).expect("default wall yields a solid");
        assert!(mesh.nodes.len() > 8, "expected layers + films + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn wall_mesh_none_for_invalid() {
        let s = InsulationWorkbenchState {
            insulation_thickness_m: 0.0,
            ..Default::default()
        };
        assert!(wall_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_insulation_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_insulation_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_insulation_workbench = true;
        run_insulation(&mut app.insulation);
        draw_workbench(&mut app);
    }
}
