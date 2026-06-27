//! The right-side **Buckling Workbench** panel — native Euler
//! column-buckling analysis over `valenx-buckling`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_buckling_workbench`,
//! toggled from the View menu. The form sets a slender column (Young's
//! modulus, smallest second moment of area, length, cross-sectional area),
//! one of the four classical end conditions, and a material yield strength;
//! "Analyze" computes the Euler critical load `P_cr = pi^2 E I / (K L)^2`,
//! the bare Euler critical stress, the J. B. Johnson parabola stress, the
//! slenderness ratio and the combined Euler/Johnson design stress and load,
//! and "Show 3-D column" loads a representative slender column solid into
//! the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_buckling::{Column, EndCondition};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Buckling Workbench.
pub struct BucklingWorkbenchState {
    /// Young's modulus `E` (Pa).
    youngs_modulus_pa: f64,
    /// Smallest second moment of area `I` (m^4).
    second_moment_area_m4: f64,
    /// Unsupported length `L` (m).
    length_m: f64,
    /// Cross-sectional area `A` (m^2).
    area_m2: f64,
    /// How the two ends are restrained (sets the effective-length factor `K`).
    end_condition: EndCondition,
    /// Material yield strength `sigma_y` (Pa), for the Johnson design branch.
    yield_strength_pa: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D column solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for BucklingWorkbenchState {
    fn default() -> Self {
        // A 3 m pinned-pinned A36-like steel strut: a 50 mm solid round bar
        // (I = 4.909e-6 m^4, A = 3.142e-3 m^2), E = 200 GPa, yield 250 MPa.
        // P_cr ~ 1.08 MN; slenderness ~ 75.9 < C_c ~ 125.7, so it is a
        // Johnson (inelastic) column and the design stress caps at ~204 MPa.
        Self {
            youngs_modulus_pa: 200.0e9,
            second_moment_area_m4: 4.909e-6,
            length_m: 3.0,
            area_m2: 3.142e-3,
            end_condition: EndCondition::PinnedPinned,
            yield_strength_pa: 250.0e6,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Buckling Workbench right-side panel. A no-op when the
/// `show_buckling_workbench` toggle is off.
pub fn draw_buckling_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_buckling_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_buckling_workbench",
        "Buckling",
        |app, ui| {
            ui.label(
                egui::RichText::new("native Euler/Johnson column buckling · valenx-buckling")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.buckling;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise).
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("Young's modulus E (Pa)");
                        ui.add(
                            egui::DragValue::new(&mut s.youngs_modulus_pa)
                                .speed(1.0e9)
                                .range(0.0..=f64::MAX),
                        )
                        .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("yield strength σy (Pa)");
                        ui.add(
                            egui::DragValue::new(&mut s.yield_strength_pa)
                                .speed(1.0e6)
                                .range(0.0..=f64::MAX),
                        )
                        .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Cross-section").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("2nd moment I (m⁴)");
                        ui.add(
                            egui::DragValue::new(&mut s.second_moment_area_m4)
                                .speed(1.0e-7)
                                .range(0.0..=f64::MAX),
                        )
                        .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("area A (m²)");
                        ui.add(
                            egui::DragValue::new(&mut s.area_m2)
                                .speed(1.0e-4)
                                .range(0.0..=f64::MAX),
                        )
                        .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("length L (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.length_m)
                                .speed(0.05)
                                .range(0.0..=f64::MAX),
                        )
                        .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("End condition").strong());
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::PinnedPinned,
                        "pinned-pinned (K = 1.0)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedPinned,
                        "fixed-pinned (K = 0.7)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedFixed,
                        "fixed-fixed (K = 0.5)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedFree,
                        "fixed-free (K = 2.0)",
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_buckling(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D column").strong())
                        .on_hover_text(
                            "Build a representative slender column (with end-restraint plates top and bottom) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Buckling").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_buckling_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.buckling` borrow is
    // released here): build the column's 3-D solid and load it.
    if app.buckling.show_3d_request {
        app.buckling.show_3d_request = false;
        load_column_3d(app);
    }
}

/// Validate the form, evaluate the column and format the readout.
fn run_buckling(s: &mut BucklingWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Column`] from the form, mapping any domain error to
/// a display string. Extracted so it is shared by the readout and the 3-D
/// gate and is unit-testable.
fn build_column(s: &BucklingWorkbenchState) -> Result<Column, String> {
    Column::new(
        s.youngs_modulus_pa,
        s.second_moment_area_m4,
        s.length_m,
        s.area_m2,
        s.end_condition,
    )
    .map_err(|e| e.to_string())
}

/// Evaluate the column and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &BucklingWorkbenchState) -> Result<String, String> {
    let col = build_column(s)?;
    let p_cr = col.critical_load();
    let sigma_cr = col.critical_stress();
    let slenderness = col.slenderness_ratio();
    let cc = col
        .transition_slenderness(s.yield_strength_pa)
        .map_err(|e| e.to_string())?;
    let johnson_stress = col
        .johnson_critical_stress(s.yield_strength_pa)
        .map_err(|e| e.to_string())?;
    let design_stress = col
        .design_critical_stress(s.yield_strength_pa)
        .map_err(|e| e.to_string())?;
    let design_load = col
        .design_critical_load(s.yield_strength_pa)
        .map_err(|e| e.to_string())?;
    let regime = if slenderness >= cc {
        "Euler (elastic)"
    } else {
        "Johnson (inelastic)"
    };

    Ok(format!(
        "end condition   : {} (K = {:.2})\n\
         E / yield σy    : {:.1} GPa / {:.1} MPa\n\
         I / area A      : {:.4e} m⁴ / {:.4e} m²\n\
         length L        : {:.3} m\n\
         effective KL    : {:.3} m\n\
         radius of gyr r : {:.3} mm\n\n\
         critical load   : {:.4} MN\n\
         Euler stress σcr: {:.2} MPa\n\
         Johnson stress  : {:.2} MPa\n\
         slenderness KL/r: {:.1}\n\
         transition C_c  : {:.1}\n\
         regime          : {}\n\
         design stress   : {:.2} MPa\n\
         design load     : {:.4} MN",
        col.end_condition.label(),
        col.factor_k(),
        s.youngs_modulus_pa / 1.0e9,
        s.yield_strength_pa / 1.0e6,
        s.second_moment_area_m4,
        s.area_m2,
        col.length,
        col.effective_length(),
        col.radius_of_gyration() * 1.0e3,
        p_cr / 1.0e6,
        sigma_cr / 1.0e6,
        johnson_stress / 1.0e6,
        slenderness,
        cc,
        regime,
        design_stress / 1.0e6,
        design_load / 1.0e6,
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

/// Build the column as a triangle [`Mesh`] — a tall thin square shaft
/// (slender in the transverse `x`/`y` directions, long in `z`) with a wider
/// end-restraint plate at the bottom and top. Representative geometry (not
/// to scale; the buckling numbers are the `valenx-buckling` result). `None`
/// for an invalid configuration.
fn column_solid_mesh(s: &BucklingWorkbenchState) -> Option<Mesh> {
    build_column(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Slender shaft (long in z, thin in x/y).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 1.5),
        Vector3::new(0.05, 0.05, 1.4),
    );
    // Bottom end-restraint plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.06),
        Vector3::new(0.3, 0.3, 0.06),
    );
    // Top end-restraint plate.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 2.94),
        Vector3::new(0.3, 0.3, 0.06),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-buckling");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D column solid and load it into the central viewport.
fn load_column_3d(app: &mut ValenxApp) {
    let Some(mesh) = column_solid_mesh(&app.buckling) else {
        app.buckling.error =
            Some("column parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<column>/valenx-buckling"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"buckling"}`** product: the representative
/// slender column (a tall thin shaft with end-restraint plates) built from the
/// canonical 3 m pinned-pinned A36-like steel strut (50 mm solid round bar, 200
/// GPa, 250 MPa yield), paired with the Euler/Johnson readout rows (critical
/// load / stress / slenderness / design stress), at a fixed 3/4 camera.
/// Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`BucklingWorkbenchState::default`].
pub(crate) fn buckling_product() -> crate::WorkspaceProduct {
    let s = BucklingWorkbenchState::default();
    let mesh = column_solid_mesh(&s).expect("canonical column ⇒ slender-column solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<column>/valenx-buckling");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical column ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Column (Euler/Johnson buckling)".into(),
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
    use std::f64::consts::PI;

    #[test]
    fn default_state_is_idle() {
        let s = BucklingWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_load_stress_and_slenderness() {
        let mut s = BucklingWorkbenchState::default();
        run_buckling(&mut s);
        assert!(
            s.error.is_none(),
            "default column should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("critical load"));
        assert!(s.result.contains("Euler stress"));
        assert!(s.result.contains("slenderness"));
        // 3 m pinned-pinned 50 mm steel round: P_cr ~ 1.0767 MN.
        assert!(s.result.contains("1.0767"));
        // Slenderness ~ 75.9 < C_c ~ 125.7, so it is a Johnson column and
        // the design stress caps below the bare Euler 342.67 MPa.
        assert!(s.result.contains("Johnson (inelastic)"));
        assert!(s.result.contains("342.67"));
    }

    #[test]
    fn analyze_rejects_zero_length() {
        let mut s = BucklingWorkbenchState {
            length_m: 0.0,
            ..Default::default()
        };
        run_buckling(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_euler_load_and_critical_stress() {
        // Ground truth: P_cr = pi^2 E I / (K L)^2 hand-computed, and the
        // critical stress is exactly P_cr / A.
        let s = BucklingWorkbenchState::default();
        let col = build_column(&s).expect("default column is valid");
        let kl = col.factor_k() * col.length; // K = 1, L = 3 -> 3.0
        let expected_pcr = PI.powi(2) * s.youngs_modulus_pa * s.second_moment_area_m4 / (kl * kl);
        assert!(
            (col.critical_load() - expected_pcr).abs() < 1.0e-3 * expected_pcr.abs(),
            "{} vs {}",
            col.critical_load(),
            expected_pcr
        );
        assert!(
            (col.critical_stress() - expected_pcr / s.area_m2).abs()
                < 1.0e-6 * (expected_pcr / s.area_m2).abs()
        );
    }

    #[test]
    fn ground_truth_johnson_parabola_stress() {
        // Ground truth: the J. B. Johnson parabola stress
        //   sigma = sy * [1 - sy (K L / r)^2 / (4 pi^2 E)]
        // hand-computed for the default 3 m pinned-pinned 50 mm steel round.
        let mut s = BucklingWorkbenchState::default();
        let col = build_column(&s).expect("default column is valid");
        let r = (s.second_moment_area_m4 / s.area_m2).sqrt(); // ~0.039527 m
        let lambda = (col.factor_k() * col.length) / r; // K=1, L=3 -> ~75.898
        let expected = s.yield_strength_pa
            * (1.0_f64
                - s.yield_strength_pa * lambda * lambda / (4.0 * PI.powi(2) * s.youngs_modulus_pa));
        // ~204.40 MPa: the crate's johnson_critical_stress matches the formula.
        let got = col
            .johnson_critical_stress(s.yield_strength_pa)
            .expect("default yield strength is positive");
        assert!(
            (got - expected).abs() < 1.0e-6 * expected.abs(),
            "Johnson stress {got} vs hand-computed {expected}"
        );
        // The readout surfaces it, formatted to 2 dp in MPa: 204.40.
        run_buckling(&mut s);
        assert!(
            s.result.contains("Johnson stress"),
            "readout should label the Johnson parabola stress: {}",
            s.result
        );
        assert!(
            s.result.contains("204.40"),
            "Johnson parabola stress should read ~204.40 MPa: {}",
            s.result
        );
    }

    #[test]
    fn fixed_fixed_carries_more_than_pinned() {
        // K = 0.5 vs 1.0 -> 4x the load (monotonic in restraint).
        let pinned = BucklingWorkbenchState::default();
        let fixed = BucklingWorkbenchState {
            end_condition: EndCondition::FixedFixed,
            ..Default::default()
        };
        let pp = build_column(&pinned).unwrap().critical_load();
        let ff = build_column(&fixed).unwrap().critical_load();
        assert!(ff > pp, "fixed-fixed should carry more: {ff} vs {pp}");
        assert!((ff / pp - 4.0).abs() < 1.0e-9);
    }

    #[test]
    fn column_mesh_for_default_is_nonempty_and_in_range() {
        let s = BucklingWorkbenchState::default();
        let mesh = column_solid_mesh(&s).expect("default column yields a solid");
        assert!(mesh.nodes.len() > 8, "expected shaft + two end plates");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn column_mesh_none_for_invalid() {
        let s = BucklingWorkbenchState {
            area_m2: 0.0,
            ..Default::default()
        };
        assert!(column_solid_mesh(&s).is_none());
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
            draw_buckling_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_buckling_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_buckling_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_buckling_workbench = true;
        run_buckling(&mut app.buckling);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The five column DragValues (E, σy, I, A, L) are SpinButtons; each
        // must be `labelled_by` its caption (egui clears a DragValue's own
        // Name), so an AI / screen reader can find the control by caption text.
        let mut app = ValenxApp::default();
        app.show_buckling_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 5,
            "expected the column numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every buckling DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["Young's modulus E (Pa)", "area A (m²)", "length L (m)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays named/invokable.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
