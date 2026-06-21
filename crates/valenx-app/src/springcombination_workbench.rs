//! The right-side **Spring Combination Workbench** panel — native
//! closed-form combinator for ideal linear (Hookean) springs over
//! `valenx-springcombination`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_springcombination_workbench`, toggled from the
//! View menu. The form takes three spring rates and a deflection, picks a
//! wiring mode (parallel or series); "Analyze" computes the equivalent
//! rate (parallel `k = sum(k_i)`, series `1/k = sum(1/k_i)`) and reports
//! the assembly force `F = k x` and stored energy `U = 0.5 k x^2`, and
//! "Show 3-D" loads representative spring solids — laid side-by-side for
//! parallel, stacked end-to-end for series — into the central viewport.

use std::f64::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_springcombination::{combine, energy_from_rate, force_from_rate, Combination, Spring};

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Spring Combination Workbench.
pub struct SpringCombinationWorkbenchState {
    /// First spring rate `k1` (N/m).
    rate1_n_per_m: f64,
    /// Second spring rate `k2` (N/m).
    rate2_n_per_m: f64,
    /// Third spring rate `k3` (N/m).
    rate3_n_per_m: f64,
    /// Assembly deflection `x` (m), used for the force / energy readout.
    deflection_m: f64,
    /// Wiring of the three springs (parallel or series).
    mode: Combination,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D spring solids (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for SpringCombinationWorkbenchState {
    fn default() -> Self {
        // Three springs of 120 / 240 / 360 N/m at a 40 mm deflection.
        // Parallel: 120+240+360 = 720 N/m, F = 28.8 N, U = 0.576 J.
        // Series:   1/k = 1/120+1/240+1/360 = 11/720 -> k = 65.45 N/m.
        Self {
            rate1_n_per_m: 120.0,
            rate2_n_per_m: 240.0,
            rate3_n_per_m: 360.0,
            deflection_m: 0.04,
            mode: Combination::Parallel,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Spring Combination Workbench right-side panel. A no-op when the
/// `show_springcombination_workbench` toggle is off.
pub fn draw_springcombination_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_springcombination_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_springcombination_workbench",
        "Spring Combination",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native closed-form linear-spring combinator · valenx-springcombination",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.springcombination;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Spring rates").strong());
                    ui.horizontal(|ui| {
                        ui.label("k₁ (N/m)");
                        ui.add(egui::DragValue::new(&mut s.rate1_n_per_m).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("k₂ (N/m)");
                        ui.add(egui::DragValue::new(&mut s.rate2_n_per_m).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("k₃ (N/m)");
                        ui.add(egui::DragValue::new(&mut s.rate3_n_per_m).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Loading").strong());
                    ui.horizontal(|ui| {
                        ui.label("deflection x (m)");
                        ui.add(egui::DragValue::new(&mut s.deflection_m).speed(0.005));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Wiring").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.mode, Combination::Parallel, "parallel");
                        ui.radio_value(&mut s.mode, Combination::Series, "series");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_springcombination(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build representative spring solids — laid side-by-side for parallel, stacked end-to-end for series — and load them into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Equivalent spring").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_springcombination_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.springcombination`
    // borrow is released here): build the spring solids and load them.
    if app.springcombination.show_3d_request {
        app.springcombination.show_3d_request = false;
        load_springs_3d(app);
    }
}

/// Validate the form, evaluate the combination and format the readout.
fn run_springcombination(s: &mut SpringCombinationWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The three validated member springs, in form order. Extracted so the
/// readout and the 3-D gate share one validation path; any non-positive or
/// non-finite rate maps to the crate's domain error.
fn springs(s: &SpringCombinationWorkbenchState) -> Result<[Spring; 3], String> {
    Ok([
        Spring::new(s.rate1_n_per_m).map_err(|e| e.to_string())?,
        Spring::new(s.rate2_n_per_m).map_err(|e| e.to_string())?,
        Spring::new(s.rate3_n_per_m).map_err(|e| e.to_string())?,
    ])
}

/// Evaluate the combination and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &SpringCombinationWorkbenchState) -> Result<String, String> {
    let members = springs(s)?;
    let k_eq = combine(s.mode, &members).map_err(|e| e.to_string())?;
    let force = force_from_rate(k_eq, s.deflection_m).map_err(|e| e.to_string())?;
    let energy = energy_from_rate(k_eq, s.deflection_m).map_err(|e| e.to_string())?;
    let mode = match s.mode {
        Combination::Parallel => "parallel",
        Combination::Series => "series",
    };

    Ok(format!(
        "member rates    : {:.1} / {:.1} / {:.1} N/m\n\
         wiring          : {mode}\n\
         deflection x    : {:.3} m\n\n\
         equivalent k    : {k_eq:.2} N/m\n\
         assembly force F: {force:.2} N\n\
         stored energy U : {energy:.3} J",
        s.rate1_n_per_m, s.rate2_n_per_m, s.rate3_n_per_m, s.deflection_m,
    ))
}

/// Append a representative coil — an axial helix swept as a thin square
/// tube — centred on `axis_origin` and running along +z, to the buffers.
/// `coils` turns of `radius`, over `height`, with a `wire` half-width.
fn push_coil(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    axis_origin: Vector3<f64>,
    radius: f64,
    height: f64,
    coils: f64,
    wire: f64,
) {
    let segments = 96_usize;
    let turns = coils * 2.0 * PI;
    // Build two rings of nodes (inner / outer offset) per station, then
    // band them into a closed quad strip along the helix.
    let mut ring_base: Vec<usize> = Vec::with_capacity(segments + 1);
    for i in 0..=segments {
        let t = i as f64 / segments as f64;
        let theta = t * turns;
        let z = axis_origin.z + t * height;
        let cx = axis_origin.x + radius * theta.cos();
        let cy = axis_origin.y + radius * theta.sin();
        // Radially-offset pair gives the tube a finite cross-section.
        let nx = theta.cos();
        let ny = theta.sin();
        ring_base.push(nodes.len());
        nodes.push(Vector3::new(cx - wire * nx, cy - wire * ny, z));
        nodes.push(Vector3::new(cx + wire * nx, cy + wire * ny, z + wire));
    }
    for i in 0..segments {
        let a = ring_base[i];
        let b = ring_base[i + 1];
        tris.extend_from_slice(&[a, a + 1, b + 1, a, b + 1, b]);
    }
}

/// Build the representative spring solids as a triangle [`Mesh`]: three
/// coils laid side-by-side (parallel — they share displacement) or stacked
/// end-to-end along the axis (series — they share force). Representative
/// geometry (not to scale; the equivalent-rate numbers are the
/// `valenx-springcombination` result). `None` for an invalid configuration.
fn springs_solid_mesh(s: &SpringCombinationWorkbenchState) -> Option<Mesh> {
    springs(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let radius = 0.18;
    let wire = 0.03;
    match s.mode {
        Combination::Parallel => {
            // Side-by-side along x, all spanning the same z extent.
            for (idx, x) in [-0.5_f64, 0.0, 0.5].into_iter().enumerate() {
                let coils = 5.0 + idx as f64; // visually distinct turn counts
                push_coil(
                    &mut nodes,
                    &mut tris,
                    Vector3::new(x, 0.0, 0.1),
                    radius,
                    0.9,
                    coils,
                    wire,
                );
            }
        }
        Combination::Series => {
            // Stacked end-to-end along z (one column).
            for (idx, z0) in [0.1_f64, 0.7, 1.3].into_iter().enumerate() {
                let coils = 4.0 + idx as f64;
                push_coil(
                    &mut nodes,
                    &mut tris,
                    Vector3::new(0.0, 0.0, z0),
                    radius,
                    0.5,
                    coils,
                    wire,
                );
            }
        }
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-springcombination");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D spring solids and load them into the central viewport.
fn load_springs_3d(app: &mut ValenxApp) {
    let Some(mesh) = springs_solid_mesh(&app.springcombination) else {
        app.springcombination.error =
            Some("spring rates are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<springs>/valenx-springcombination"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical springcombination workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn springcombination_product() -> crate::WorkspaceProduct {
    let s = SpringCombinationWorkbenchState::default();
    let mesh =
        springs_solid_mesh(&s).expect("canonical springcombination ⇒ spring-set solid builds");
    let loaded =
        crate::products_registry::loaded_mesh_from(mesh, "<springcombination>/valenx-springs");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical springcombination ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Spring combination (series/parallel)".into(),
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
        let s = SpringCombinationWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_parallel_reports_equivalent_force_and_energy() {
        let mut s = SpringCombinationWorkbenchState::default();
        run_springcombination(&mut s);
        assert!(
            s.error.is_none(),
            "default springs should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("equivalent k"));
        assert!(s.result.contains("assembly force F"));
        assert!(s.result.contains("stored energy U"));
        // Parallel 120+240+360 = 720 N/m; F = 720*0.04 = 28.8 N;
        // U = 0.5*720*0.04^2 = 0.576 J.
        assert!(s.result.contains("720.00"));
        assert!(s.result.contains("28.80"));
        assert!(s.result.contains("0.576"));
    }

    #[test]
    fn analyze_series_softens_below_smallest_member() {
        let mut s = SpringCombinationWorkbenchState {
            mode: Combination::Series,
            ..Default::default()
        };
        run_springcombination(&mut s);
        assert!(s.error.is_none(), "series should analyze: {:?}", s.error);
        // 1/k = 1/120 + 1/240 + 1/360 = 11/720 -> k = 65.4545... N/m.
        assert!(s.result.contains("65.45"));
    }

    #[test]
    fn analyze_rejects_non_positive_rate() {
        let mut s = SpringCombinationWorkbenchState {
            rate2_n_per_m: 0.0,
            ..Default::default()
        };
        run_springcombination(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn two_equal_springs_ground_truth_parallel_2k_series_k_over_2() {
        // Ground truth (hand-computed): two equal springs of rate k give
        // 2k in parallel and k/2 in series. Use a degenerate third member
        // that is the identity for each rule: a huge rate is ~open in
        // series (1/k -> 0) and we compare the pure two-spring rules
        // directly through the crate combinator.
        let k = 320.0;
        let pair = [Spring::new(k).unwrap(), Spring::new(k).unwrap()];
        let parallel = combine(Combination::Parallel, &pair).unwrap();
        let series = combine(Combination::Series, &pair).unwrap();
        // Pin the exact literal before .abs() to keep the type as f64.
        assert!((parallel - 2.0 * k).abs() < 1.0e-9_f64);
        assert!((series - k / 2.0).abs() < 1.0e-9_f64);
    }

    #[test]
    fn springs_mesh_for_default_is_nonempty_and_in_range() {
        for mode in [Combination::Parallel, Combination::Series] {
            let s = SpringCombinationWorkbenchState {
                mode,
                ..Default::default()
            };
            let mesh = springs_solid_mesh(&s).expect("default springs yield a solid");
            assert!(mesh.nodes.len() > 8, "expected three swept coils");
            let n = mesh.nodes.len() as u32;
            for blk in &mesh.element_blocks {
                assert!(!blk.connectivity.is_empty());
                assert_eq!(blk.connectivity.len() % 3, 0);
                assert!(blk.connectivity.iter().all(|&i| i < n));
            }
        }
    }

    #[test]
    fn springs_mesh_none_for_invalid() {
        let s = SpringCombinationWorkbenchState {
            rate1_n_per_m: -1.0,
            ..Default::default()
        };
        assert!(springs_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_springcombination_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_springcombination_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_springcombination_workbench = true;
        run_springcombination(&mut app.springcombination);
        draw_workbench(&mut app);
    }
}
