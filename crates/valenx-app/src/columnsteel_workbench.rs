//! The right-side **Steel Column Workbench** panel — native Euler-Johnson
//! axial-compression buckling over `valenx-columnsteel`.
//!
//! Mirrors the Beam / Heat Transfer workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_columnsteel_workbench`,
//! toggled from the View menu. The form drives a
//! [`valenx_columnsteel::Column`] built from geometry (Young's modulus,
//! yield stress, effective-length factor `K` from the chosen end condition,
//! unbraced length and radius of gyration). "Analyze" reports the
//! slenderness ratio `λ = K L / r`, the column-slenderness transition
//! `Cc`, the governing regime (elastic Euler vs inelastic Johnson), the
//! elastic Euler reference stress `Fe = π² E / λ²`, the critical buckling
//! stress `Fcr`, and the AISC-ASD allowable stress and axial load;
//! "Show 3-D column" loads a representative I-section column solid into
//! the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_columnsteel::Column;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Idealized end-condition support, selecting the standard AISC theoretical
/// effective-length factor `K` fed into [`Column::from_geometry`].
///
/// The crate itself does not derive `K` from end conditions — its docs are
/// explicit that the caller supplies `K` — so this enum is purely the
/// workbench's convenience mapping to the textbook theoretical values
/// (AISC Commentary Table C-A-7.1).
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum EndCondition {
    /// Both ends pinned (theoretical `K = 1.0`).
    PinnedPinned,
    /// Both ends fixed (theoretical `K = 0.5`).
    FixedFixed,
    /// One end fixed, one end pinned (theoretical `K = 0.7`).
    FixedPinned,
    /// One end fixed, one end free — a flagpole (theoretical `K = 2.0`).
    FixedFree,
}

impl EndCondition {
    /// The theoretical effective-length factor `K` for this end condition.
    fn k(self) -> f64 {
        match self {
            EndCondition::PinnedPinned => 1.0,
            EndCondition::FixedFixed => 0.5,
            EndCondition::FixedPinned => 0.7,
            EndCondition::FixedFree => 2.0,
        }
    }
}

/// Persistent form + result state for the Steel Column Workbench.
pub struct ColumnSteelWorkbenchState {
    /// Young's modulus `E` (Pa).
    youngs_modulus: f64,
    /// Yield stress `Fy` (Pa).
    yield_stress: f64,
    /// Idealized end condition, selecting the effective-length factor `K`.
    end_condition: EndCondition,
    /// Unbraced length `L` (m).
    length_m: f64,
    /// Least radius of gyration `r` (m).
    radius_gyration_m: f64,
    /// Gross cross-sectional area `A` (m^2).
    area_m2: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D column solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ColumnSteelWorkbenchState {
    fn default() -> Self {
        // A992 W-shape steel column in SI: E = 200 GPa, Fy = 345 MPa,
        // pinned-pinned (K = 1.0), 4.0 m unbraced, r = 0.060 m,
        // A = 0.013 m^2. lambda = 66.667 < Cc = 106.972 -> Johnson;
        // Fcr ~ 278 MPa, Fa ~ 149 MPa, P_allow ~ 1932 kN.
        Self {
            youngs_modulus: 200.0e9,
            yield_stress: 345.0e6,
            end_condition: EndCondition::PinnedPinned,
            length_m: 4.0,
            radius_gyration_m: 0.06,
            area_m2: 0.013,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Steel Column Workbench right-side panel. A no-op when the
/// `show_columnsteel_workbench` toggle is off.
pub fn draw_columnsteel_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_columnsteel_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_columnsteel_workbench",
        "Steel Column",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Euler-Johnson (AISC-ASD) column buckling · valenx-columnsteel",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.columnsteel;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Material").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise).
                    ui.horizontal(|ui| {
                        let l = ui.label("Young's E (Pa)");
                        ui.add(egui::DragValue::new(&mut s.youngs_modulus).speed(1.0e9))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("yield Fy (Pa)");
                        ui.add(egui::DragValue::new(&mut s.yield_stress).speed(1.0e6))
                            .labelled_by(l.id);
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("End condition (K)").strong());
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::PinnedPinned,
                        "pinned-pinned (K=1.0)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedFixed,
                        "fixed-fixed (K=0.5)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedPinned,
                        "fixed-pinned (K=0.7)",
                    );
                    ui.radio_value(
                        &mut s.end_condition,
                        EndCondition::FixedFree,
                        "fixed-free / flagpole (K=2.0)",
                    );

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    ui.horizontal(|ui| {
                        let l = ui.label("length L (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(0.05))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("radius of gyration r (m)");
                        ui.add(egui::DragValue::new(&mut s.radius_gyration_m).speed(0.002))
                            .labelled_by(l.id);
                    });
                    ui.horizontal(|ui| {
                        let l = ui.label("area A (m²)");
                        ui.add(egui::DragValue::new(&mut s.area_m2).speed(0.001))
                            .labelled_by(l.id);
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_columnsteel(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D column").strong())
                        .on_hover_text(
                            "Build a representative I-section steel column (two flanges + a web) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Capacity").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_columnsteel_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.columnsteel` borrow is
    // released here): build the column's 3-D solid and load it.
    if app.columnsteel.show_3d_request {
        app.columnsteel.show_3d_request = false;
        load_column_3d(app);
    }
}

/// Validate the form, evaluate the column and format the readout.
fn run_columnsteel(s: &mut ColumnSteelWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Column`] for the current form, the gate shared by
/// the readout and the 3-D solid. Extracted so it is unit-testable.
fn column_of(s: &ColumnSteelWorkbenchState) -> Result<Column, String> {
    Column::from_geometry(
        s.youngs_modulus,
        s.yield_stress,
        s.end_condition.k(),
        s.length_m,
        s.radius_gyration_m,
    )
    .map_err(|e| e.to_string())
}

/// Evaluate the column and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ColumnSteelWorkbenchState) -> Result<String, String> {
    let col = column_of(s)?;
    let lambda = col.slenderness();
    let cc = col.cc();
    let regime = col.regime().as_str();
    let fcr = col.critical_stress();
    // Elastic Euler reference stress Fe = π² E / λ², the bare elastic
    // buckling curve. For a long (Euler-regime) column this equals the
    // governing Fcr; for a short / intermediate (Johnson-regime) column it
    // is the higher elastic prediction the inelastic Johnson parabola
    // corrects downward, so showing both makes the regime tangible.
    let fe = col.euler_stress_unchecked();
    let fs = col.factor_of_safety_aisc();
    let fa = col.allowable_stress();
    let p_allow = col.allowable_load(s.area_m2).map_err(|e| e.to_string())?;

    Ok(format!(
        "E / Fy          : {:.0} / {:.0} MPa\n\
         end condition K : {:.2}\n\
         length / r      : {:.3} m / {:.4} m\n\
         area A          : {:.5} m²\n\n\
         slenderness λ   : {:.3}\n\
         transition Cc   : {:.3}\n\
         regime          : {regime}\n\
         Euler ref. Fe   : {:.2} MPa\n\
         critical Fcr    : {:.2} MPa\n\
         safety factor   : {:.3}\n\
         allowable Fa    : {:.2} MPa\n\
         allowable load  : {:.1} kN",
        s.youngs_modulus / 1.0e6,
        s.yield_stress / 1.0e6,
        s.end_condition.k(),
        s.length_m,
        s.radius_gyration_m,
        s.area_m2,
        lambda,
        cc,
        fe / 1.0e6,
        fcr / 1.0e6,
        fs,
        fa / 1.0e6,
        p_allow / 1.0e3,
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

/// Build the steel column as a triangle [`Mesh`] — a tall I-section
/// (wide-flange) member approximated by three vertical boxes: two flanges
/// and a connecting web, on a small base pad. Representative geometry (not
/// to scale; the capacity numbers are the `valenx-columnsteel` result).
/// `None` for an invalid configuration.
fn column_solid_mesh(s: &ColumnSteelWorkbenchState) -> Option<Mesh> {
    column_of(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Tall I-section: half-height 1.6 (z runs 0..3.2), flange half-width
    // 0.25, web half-thickness 0.04, flange half-thickness 0.06.
    let half_h = 1.6;
    let zc = half_h + 0.08; // sit on top of the base pad
    let flange_off = 0.31; // ±x position of the two flanges

    // Web (the central plate).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, zc),
        Vector3::new(0.04, 0.25, half_h),
    );
    // Top flange (+x).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(flange_off, 0.0, zc),
        Vector3::new(0.06, 0.25, half_h),
    );
    // Bottom flange (-x).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-flange_off, 0.0, zc),
        Vector3::new(0.06, 0.25, half_h),
    );
    // Base pad.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.04),
        Vector3::new(0.5, 0.4, 0.04),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-columnsteel");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D column solid and load it into the central viewport.
fn load_column_3d(app: &mut ValenxApp) {
    let Some(mesh) = column_solid_mesh(&app.columnsteel) else {
        app.columnsteel.error =
            Some("column parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<column>/valenx-columnsteel"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical columnsteel workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn columnsteel_product() -> crate::WorkspaceProduct {
    let s = ColumnSteelWorkbenchState::default();
    let mesh = column_solid_mesh(&s).expect("canonical columnsteel ⇒ wide-flange solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<columnsteel>/valenx-column");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical columnsteel ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Steel column (AISC compression)".into(),
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
        let s = ColumnSteelWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_slenderness_regime_and_capacity() {
        let mut s = ColumnSteelWorkbenchState::default();
        run_columnsteel(&mut s);
        assert!(
            s.error.is_none(),
            "default column should analyze: {:?}",
            s.error
        );
        // lambda = 1.0 * 4.0 / 0.06 = 66.667, Cc = 106.972, Johnson.
        assert!(s.result.contains("66.667"));
        assert!(s.result.contains("106.972"));
        assert!(s.result.contains("johnson"));
        // Fcr ~ 278.00 MPa, Fa ~ 148.65 MPa.
        assert!(s.result.contains("278.00"));
        assert!(s.result.contains("148.65"));
    }

    #[test]
    fn analyze_rejects_zero_radius_of_gyration() {
        let mut s = ColumnSteelWorkbenchState {
            radius_gyration_m: 0.0,
            ..Default::default()
        };
        run_columnsteel(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn slenderness_ground_truth_k_l_over_r_and_euler_for_slender() {
        // Ground truth 1: slenderness lambda = K * L / r, computed by hand
        // for the pinned-pinned default geometry.
        let s = ColumnSteelWorkbenchState::default();
        let col = column_of(&s).expect("default column is valid");
        let lambda_hand: f64 = 1.0 * 4.0 / 0.06;
        assert!((col.slenderness() - lambda_hand).abs() < 1e-9);

        // Ground truth 2: for a deliberately slender column (lambda well
        // past Cc) the critical stress is the elastic Euler curve
        // Fcr = pi^2 * E / lambda^2, hand-computed.
        let e: f64 = 200.0e9;
        let fy: f64 = 345.0e6;
        let lambda: f64 = 150.0;
        let slender = Column::new(e, fy, lambda).expect("slender column is valid");
        assert_eq!(slender.regime().as_str(), "euler");
        let fcr_hand = PI * PI * e / (lambda * lambda);
        assert!((slender.critical_stress() - fcr_hand).abs() < 1e-3);
    }

    #[test]
    fn euler_reference_stress_ground_truth_and_in_readout() {
        // The readout surfaces the elastic Euler reference stress
        // Fe = pi^2 * E / lambda^2. For the Johnson-regime default it is
        // the higher elastic prediction that the inelastic Johnson curve
        // corrects downward, so it must differ from the governing Fcr.
        let mut s = ColumnSteelWorkbenchState::default();
        run_columnsteel(&mut s);
        assert!(
            s.error.is_none(),
            "default column should analyze: {:?}",
            s.error
        );

        // Ground truth: lambda = 1.0 * 4.0 / 0.06 = 66.6667;
        // Fe = pi^2 * 200e9 / 66.6667^2 = 4.44132e8 Pa = 444.13 MPa.
        let e: f64 = 200.0e9;
        let lambda: f64 = 1.0 * 4.0 / 0.06;
        let fe_hand_mpa = PI * PI * e / (lambda * lambda) / 1.0e6;
        let one = 1.0_f64;
        assert!(
            (fe_hand_mpa - 444.132_198).abs() < one * 1e-3,
            "hand Fe = {fe_hand_mpa} MPa"
        );

        // The crate accessor agrees with the hand value.
        let col = column_of(&s).expect("default column is valid");
        assert!((col.euler_stress_unchecked() / 1.0e6 - fe_hand_mpa).abs() < 1e-6);

        // The formatted readout carries the Euler reference line, and it is
        // a genuinely distinct value from the Johnson Fcr (278.00 MPa).
        assert!(s.result.contains("Euler ref. Fe"));
        assert!(s.result.contains("444.13"));
        assert!(s.result.contains("278.00"));
    }

    #[test]
    fn end_condition_k_scales_slenderness() {
        // Fixed-fixed (K = 0.5) halves the slenderness vs pinned-pinned
        // (K = 1.0) for the same geometry.
        let pinned = ColumnSteelWorkbenchState::default();
        let fixed = ColumnSteelWorkbenchState {
            end_condition: EndCondition::FixedFixed,
            ..Default::default()
        };
        let lp = column_of(&pinned).unwrap().slenderness();
        let lf = column_of(&fixed).unwrap().slenderness();
        assert!((lf - 0.5 * lp).abs() < 1e-9, "K=0.5 should halve lambda");
    }

    #[test]
    fn column_mesh_for_default_is_nonempty_and_in_range() {
        let s = ColumnSteelWorkbenchState::default();
        let mesh = column_solid_mesh(&s).expect("default column yields a solid");
        assert!(mesh.nodes.len() > 8, "expected web + 2 flanges + base");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn column_mesh_none_for_invalid() {
        let s = ColumnSteelWorkbenchState {
            radius_gyration_m: 0.0,
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

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_columnsteel_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_columnsteel_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_columnsteel_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_columnsteel_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_columnsteel_workbench = true;
        run_columnsteel(&mut app.columnsteel);
        app.columnsteel.error = Some("invalid column parameters".to_string());
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The material + geometry DragValues are SpinButtons; each must be
        // `labelled_by` its caption (egui clears a DragValue's own Name), so an
        // AI / screen reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_columnsteel_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // E, Fy, length, radius of gyration, area.
        assert!(
            spin_buttons.len() >= 5,
            "expected the column numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every column DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["Young's E (Pa)", "length L (m)", "area A (m²)"] {
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
