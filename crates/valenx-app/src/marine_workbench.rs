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

use valenx_marine::{Hull, FRESHWATER_DENSITY, SEAWATER_DENSITY};
use valenx_mesh::Mesh;

use crate::mesh_prims::MeshBuilder;
use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Antifouling red for the hull below the waterline.
const ANTIFOUL: [f32; 3] = [0.55, 0.13, 0.13];
/// Topside grey for the hull above the waterline.
const TOPSIDE: [f32; 3] = [0.40, 0.43, 0.48];
/// Teak-ish deck.
const DECK: [f32; 3] = [0.62, 0.50, 0.34];
/// Dark keel fin.
const KEEL: [f32; 3] = [0.18, 0.18, 0.20];

/// Number of transverse stations lofted along the hull length (stern → bow).
const STATIONS: usize = 15;
/// Half-section points from keel/waterline up one side (per band).
const SECTION_POINTS: usize = 8;

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

/// The longitudinal **fullness** of the hull at fractional length `u ∈ [0, 1]`
/// (0 = stern, 1 = bow). 1.0 amidships, tapering toward both ends; the bow
/// (forward) is finer than the stern. Raised to an exponent driven by the block
/// coefficient `cb` so a fuller hull (high Cb) keeps its sections fuller for
/// longer, a fine hull (low Cb) pinches in earlier — keeping the rendered form
/// consistent with the hydrostatics readout's `Cb`.
fn fullness(u: f64, cb: f64) -> f64 {
    // Base bell shape, peak at amidships (u≈0.5), zero at the very ends.
    let bell = (std::f64::consts::PI * u).sin();
    // Bias the peak slightly aft so the bow is finer than the stern.
    let aft_bias = 1.0 - 0.18 * (u - 0.5);
    // Sharpness: fuller hull ⇒ flatter top (smaller exponent).
    let sharp = (1.6 - cb).clamp(0.5, 1.4);
    (bell.powf(sharp) * aft_bias).clamp(0.0, 1.0)
}

/// One transverse **half-section ring band** at fractional length `u`, spanning
/// the vertical band from `z_lo` to `z_hi` (a fraction of the local depth). The
/// section is a closed loop: down the **port** side from the upper edge to the
/// keel/lower edge, across to **starboard**, up the starboard side, then closed
/// back across the top — so two adjacent bands loft into a watertight hull
/// skin. `hl`/`hb` are half-length/half-beam, `depth` the keel-to-waterline
/// draft, `f` the local [`fullness`]. The section narrows in beam and rises off
/// the keel toward the ends via `f`, giving the fine entrance / run.
fn hull_band(
    u: f64,
    z_lo: f64,
    z_hi: f64,
    hl: f64,
    hb: f64,
    depth: f64,
    f: f64,
) -> Vec<[f64; 3]> {
    // Longitudinal x: stern (−hl) at u=0 to a raked bow (+hl) at u=1.
    let x = -hl + 2.0 * hl * u;
    // Local half-beam at this station (full amidships, → ~0 at the ends).
    let local_hb = hb * (0.12 + 0.88 * f);
    // Keel rise: the bottom lifts toward the ends (rocker / forefoot).
    let keel_z = (1.0 - f) * depth * 0.55;
    let n = SECTION_POINTS;
    // Half-beam at vertical fraction `t` (0 top → 1 bottom): rounds in toward
    // the keel (a gentle bilge curve), never quite zero.
    let half_y = |t: f64| -> f64 { local_hb * (1.0 - t * t * 0.85).max(0.06) };
    let z_at = |t: f64| -> f64 { (z_hi + (z_lo - z_hi) * t).max(keel_z) };
    let mut ring = Vec::with_capacity(2 * n);
    // Port side: upper edge (t=0) → lower edge (t=1).
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        ring.push([x, -half_y(t), z_at(t)]);
    }
    // Starboard side: lower edge → upper edge (mirror), skipping the shared
    // bottom keel point (i=n-1) so the closed loop has no duplicate vertex but
    // keeps both top corners at full beam.
    for i in (0..n - 1).rev() {
        let t = i as f64 / (n - 1) as f64;
        ring.push([x, half_y(t), z_at(t)]);
    }
    ring
}

/// Build a representative ship hull as a triangle [`Mesh`] **with per-vertex
/// colours** — a **lofted hull** skinned from [`STATIONS`] transverse station
/// sections along the length (fuller amidships, fine at the bow/stern entrance
/// and run, with rocker lifting the keel toward the ends), split into a
/// **below-waterline** band (antifoul red, keel → waterline) and a **topside**
/// band (grey, waterline → sheer), plus a flat **deck** and a centreline
/// **keel** fin. The station beam/draft are driven by the workbench's own
/// length/beam/draft and a [`fullness`] curve shaped by `Cb`, so the rendered
/// form tracks the hydrostatics readout. Length runs along x, beam along y,
/// draft from the keel `z = 0` up to the waterline `z = T`. The reported
/// hydrostatics still use the box-form `Cb` model from `valenx-marine`. `None`
/// for an invalid hull.
///
/// Returns `(mesh, colors)` with `colors.len() == 3 × triangle_count`, ready
/// for [`crate::WorkspaceProduct::vertex_colors`].
fn hull_solid_mesh_colored(s: &MarineWorkbenchState) -> Option<(Mesh, Vec<[f32; 3]>)> {
    let hull = build_hull(s).ok()?;
    let (hl, hb, t) = (hull.length_m / 2.0, hull.beam_m / 2.0, hull.draft_m);
    let cb = hull.block_coefficient;
    // Freeboard: topside rises above the waterline by ~40 % of the draft.
    let freeboard = t * 0.4;
    let deck_z = t + freeboard;

    // Build the station sections for the two vertical bands. Endpoints (u=0,
    // u=1) get a tiny non-zero fullness so the stem/transom rings stay valid
    // closed loops (capped) rather than collapsing to a line.
    let mut below: Vec<Vec<[f64; 3]>> = Vec::with_capacity(STATIONS);
    let mut top: Vec<Vec<[f64; 3]>> = Vec::with_capacity(STATIONS);
    for k in 0..STATIONS {
        let u = k as f64 / (STATIONS - 1) as f64;
        let f = fullness(u, cb).max(0.08);
        below.push(hull_band(u, 0.0, t, hl, hb, t, f));
        top.push(hull_band(u, t, deck_z, hl, hb, t, f));
    }

    let mut b = MeshBuilder::new();
    // Below-waterline skin (red) and topside skin (grey), each capped so the
    // transom (stern) and the bow stem are closed.
    b.loft(&below, true, ANTIFOUL);
    b.loft(&top, true, TOPSIDE);

    // Deck: a flat slab spanning the full length, slightly inset, at the sheer.
    let deck_hb = hb * 0.96;
    b.cuboid(
        [0.0, 0.0, deck_z],
        [2.0 * hl * 0.98, 2.0 * deck_hb, t * 0.08],
        DECK,
    );

    // Keel: a thin centreline fin running most of the length below the hull.
    b.cuboid(
        [hl * -0.05, 0.0, -t * 0.12],
        [2.0 * hl * 0.6, hb * 0.06, t * 0.24],
        KEEL,
    );

    let (mut mesh, colors) = b.into_mesh_and_colors();
    mesh.id = "valenx-marine-hull".to_string();
    Some((mesh, colors))
}

/// Build the hull [`Mesh`] (without the colour metadata) for the central
/// viewport. See [`hull_solid_mesh_colored`].
fn hull_solid_mesh(s: &MarineWorkbenchState) -> Option<Mesh> {
    hull_solid_mesh_colored(s).map(|(mesh, _colors)| mesh)
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

/// The agent-bridge **`show_3d{kind:"marine"}`** product: the canonical ship
/// hull built as a 3-D solid, paired with the workbench's own hydrostatics
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`MarineWorkbenchState::default`]. This workbench
/// formats its readout via [`run_marine`] into `s.result` (it has no separate
/// `compute()`), so the builder runs the analysis first and reads that field.
pub(crate) fn marine_product() -> crate::WorkspaceProduct {
    let mut s = MarineWorkbenchState::default();
    let (mesh, colors) = hull_solid_mesh_colored(&s).expect("canonical marine ⇒ hull solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<marine>/valenx-hull");
    run_marine(&mut s);
    let lines = crate::products_registry::lines_from_readout(&s.result);
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Marine hull (hydrostatics)".into(),
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
        // The lofted hull has far more than the old 10-node block (two skinned
        // bands of STATIONS sections + deck + keel).
        assert!(
            mesh.nodes.len() > 100,
            "lofted hull has many station vertices, got {}",
            mesh.nodes.len()
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
        // Length runs along x (≈ full L), beam along y (≈ full B), keel/deck
        // span the draft + freeboard in z.
        let max_x = mesh.nodes.iter().map(|p| p.x).fold(f64::MIN, f64::max);
        let max_y = mesh.nodes.iter().map(|p| p.y).fold(f64::MIN, f64::max);
        assert!((max_x - s.length_m / 2.0).abs() < s.length_m * 0.05, "spans L");
        assert!(max_y <= s.beam_m / 2.0 + 1e-6, "within the beam");
    }

    #[test]
    fn hull_mesh_none_for_invalid_hull() {
        let s = MarineWorkbenchState {
            draft_m: 0.0,
            ..Default::default()
        };
        assert!(hull_solid_mesh(&s).is_none());
    }

    #[test]
    fn fullness_is_peaked_amidships_and_fine_at_ends() {
        // Full amidships, fine at the bow/stern; the bow (u→1) is finer than the
        // stern (u→0) at symmetric offsets.
        let cb = 0.7;
        assert!(fullness(0.5, cb) > 0.85, "fullest amidships");
        assert!(fullness(0.02, cb) < 0.3, "fine at the stern");
        assert!(fullness(0.98, cb) < 0.3, "fine at the bow");
        assert!(
            fullness(0.75, cb) < fullness(0.25, cb),
            "bow finer than stern (aft bias)"
        );
        // A fuller hull (higher Cb) stays fuller off-amidships than a fine one.
        assert!(fullness(0.3, 0.85) > fullness(0.3, 0.5), "high Cb keeps fullness");
    }

    #[test]
    fn hull_band_is_a_closed_within_envelope_ring() {
        // A station ring is a closed loop of 2·n−1 points (port n + starboard
        // n−1, sharing the keel point), symmetric in y, with x at the station
        // and z within the [keel, deck] band.
        let ring = hull_band(0.5, 0.0, 6.0, 60.0, 10.0, 6.0, 1.0);
        assert_eq!(ring.len(), 2 * SECTION_POINTS - 1);
        let max_y = ring.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
        let min_y = ring.iter().map(|p| p[1]).fold(f64::MAX, f64::min);
        assert!((max_y + min_y).abs() < 1e-9, "port/starboard symmetric");
        assert!(max_y > 0.0 && max_y <= 10.0 + 1e-9, "within the half-beam");
        let min_z = ring.iter().map(|p| p[2]).fold(f64::MAX, f64::min);
        assert!(min_z >= -1e-9, "keel at or above z=0");
    }

    #[test]
    fn hull_carries_vertex_aligned_colours() {
        // The two lofted bands + deck + keel ship per-vertex colours aligned to
        // the renderer's coloured path (3 / triangle), with the antifoul,
        // topside and deck colours all present.
        let s = MarineWorkbenchState::default();
        let (mesh, colors) = hull_solid_mesh_colored(&s).expect("default hull builds coloured");
        assert!(!mesh.nodes.is_empty(), "non-empty mesh");
        assert!(mesh.total_elements() > 0, "mesh has triangles");
        assert_eq!(
            colors.len(),
            mesh.total_elements() * 3,
            "vertex_colors must equal 3 × triangle count"
        );
        assert!(colors.contains(&ANTIFOUL), "below-waterline colour present");
        assert!(colors.contains(&TOPSIDE), "topside colour present");
        assert!(colors.contains(&DECK), "deck colour present");
        for c in &colors {
            for ch in c {
                assert!(ch.is_finite() && (0.0..=1.0).contains(ch));
            }
        }
    }

    #[test]
    fn hull_product_is_coloured_and_aligned() {
        let product = marine_product();
        let loaded = product.mesh.as_ref().expect("marine product has a mesh");
        let colors = product
            .vertex_colors
            .as_ref()
            .expect("marine product carries vertex_colors");
        assert_eq!(
            colors.len(),
            loaded.mesh.total_elements() * 3,
            "product colours aligned to the coloured path"
        );
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
