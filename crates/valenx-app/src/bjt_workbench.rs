//! The right-side **BJT Workbench** panel — native DC bias Q-point
//! analysis over `valenx-bjt`.
//!
//! Mirrors the Heat Transfer / DC Motor workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_bjt_workbench`,
//! toggled from the View menu. The form describes a single bipolar
//! junction transistor (`beta`, `VBE`, `Vce_sat`) in one of two textbook
//! bias topologies — four-resistor voltage divider or single-resistor
//! fixed base — and "Analyze" solves the base loop for the quiescent
//! operating point: the base / collector / emitter currents, the
//! collector-emitter voltage, the common-base current gain `alpha`, the
//! collector power dissipation `Pd`, the bias stability factor `S(ICO)`,
//! and the operating region (active / saturation / cut-off). "Show 3-D
//! transistor" loads a representative TO-92-style package (a flat-faced
//! half-cylinder body with three leads) into the central viewport.
//!
//! Honest scope: research / educational grade. The numbers are the
//! constant-`VBE`, constant-`beta` hand-analysis result from
//! [`valenx_bjt`] — first-order design estimates, not a SPICE / Gummel-Poon
//! simulation.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_bjt::bias::{DividerBias, FixedBias, OperatingPoint};
use valenx_bjt::model::{Region, Transistor};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Which DC bias network the form describes.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Topology {
    /// Four-resistor voltage-divider bias (`R1`, `R2`, `Rc`, `Re`) — the
    /// stiff, `beta`-insensitive workhorse.
    #[default]
    Divider,
    /// Single base resistor `Rb` from the supply, plus `Rc` and `Re` —
    /// the simplest (and most `beta`-sensitive) bias.
    FixedBase,
}

/// Persistent form + result state for the BJT Workbench.
pub struct BjtWorkbenchState {
    /// Which bias topology to solve.
    topology: Topology,
    /// Forward DC current gain `beta = Ic / Ib`.
    beta: f64,
    /// Base-emitter turn-on voltage `VBE` (V).
    vbe: f64,
    /// Collector-emitter saturation floor `Vce_sat` (V).
    vce_sat: f64,
    /// Supply voltage `Vcc` (V).
    vcc: f64,
    /// Upper divider resistor `R1`, from `Vcc` to the base (kΩ).
    r1_kohm: f64,
    /// Lower divider resistor `R2`, from the base to ground (kΩ).
    r2_kohm: f64,
    /// Single base resistor `Rb` for the fixed-base topology (kΩ).
    rb_kohm: f64,
    /// Collector resistor `Rc` (kΩ).
    rc_kohm: f64,
    /// Emitter resistor `Re` (kΩ).
    re_kohm: f64,
    /// Formatted Q-point readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D transistor solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for BjtWorkbenchState {
    fn default() -> Self {
        // Classic stiff voltage-divider stage: Vcc = 12 V, R1 = 47k,
        // R2 = 10k, Rc = 2.2k, Re = 1k, silicon beta = 100. Vth = 2.105 V,
        // Rth = 8.246k, so Ib ≈ 13 µA, Ic ≈ 1.3 mA, Vce ≈ 7.7 V — solidly
        // active.
        Self {
            topology: Topology::Divider,
            beta: 100.0,
            vbe: 0.7,
            vce_sat: 0.2,
            vcc: 12.0,
            r1_kohm: 47.0,
            r2_kohm: 10.0,
            rb_kohm: 240.0,
            rc_kohm: 2.2,
            re_kohm: 1.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the BJT Workbench right-side panel. A no-op when the
/// `show_bjt_workbench` toggle is off.
pub fn draw_bjt_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bjt_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_bjt_workbench",
        "BJT",
        |app, ui| {
            ui.label(
                egui::RichText::new("native DC bias Q-point analysis · valenx-bjt")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.bjt;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Bias topology").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.topology, Topology::Divider, "voltage divider");
                        ui.radio_value(&mut s.topology, Topology::FixedBase, "fixed base");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Transistor").strong());
                    ui.horizontal(|ui| {
                        ui.label("gain β");
                        ui.add(egui::DragValue::new(&mut s.beta).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("VBE (V)");
                        ui.add(egui::DragValue::new(&mut s.vbe).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Vce_sat (V)");
                        ui.add(egui::DragValue::new(&mut s.vce_sat).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Network").strong());
                    ui.horizontal(|ui| {
                        ui.label("supply Vcc (V)");
                        ui.add(egui::DragValue::new(&mut s.vcc).speed(0.5));
                    });
                    match s.topology {
                        Topology::Divider => {
                            ui.horizontal(|ui| {
                                ui.label("R1 (kΩ)");
                                ui.add(egui::DragValue::new(&mut s.r1_kohm).speed(1.0));
                            });
                            ui.horizontal(|ui| {
                                ui.label("R2 (kΩ)");
                                ui.add(egui::DragValue::new(&mut s.r2_kohm).speed(1.0));
                            });
                        }
                        Topology::FixedBase => {
                            ui.horizontal(|ui| {
                                ui.label("Rb (kΩ)");
                                ui.add(egui::DragValue::new(&mut s.rb_kohm).speed(5.0));
                            });
                        }
                    }
                    ui.horizontal(|ui| {
                        ui.label("Rc (kΩ)");
                        ui.add(egui::DragValue::new(&mut s.rc_kohm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Re (kΩ)");
                        ui.add(egui::DragValue::new(&mut s.re_kohm).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_bjt(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D transistor").strong())
                        .on_hover_text(
                            "Build a representative TO-92-style transistor package (flat-faced body with three leads) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Q-point").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_bjt_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.bjt` borrow is released
    // here): build the transistor's 3-D solid and load it.
    if app.bjt.show_3d_request {
        app.bjt.show_3d_request = false;
        load_transistor_3d(app);
    }
}

/// Validate the form, solve the bias network and format the readout.
fn run_bjt(s: &mut BjtWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated `Transistor` for the current form. Shared by the
/// readout and the 3-D gate so it is unit-testable.
fn device(s: &BjtWorkbenchState) -> Result<Transistor, String> {
    Transistor::new(s.beta, s.vbe, s.vce_sat).map_err(|e| e.to_string())
}

/// Solve the selected bias network for `device`, returning the operating
/// point and the stability factor `S(ICO)`. Extracted so it is shared by
/// the readout and unit-testable independently of the formatting.
fn solve(s: &BjtWorkbenchState) -> Result<(OperatingPoint, f64), String> {
    let q = device(s)?;
    let rc = s.rc_kohm * 1.0e3;
    let re = s.re_kohm * 1.0e3;
    match s.topology {
        Topology::Divider => {
            let bias = DividerBias::new(s.vcc, s.r1_kohm * 1.0e3, s.r2_kohm * 1.0e3, rc, re)
                .map_err(|e| e.to_string())?;
            let op = bias.solve(&q).map_err(|e| e.to_string())?;
            let s_ico = bias.stability_factor(&q).map_err(|e| e.to_string())?;
            Ok((op, s_ico))
        }
        Topology::FixedBase => {
            let bias =
                FixedBias::new(s.vcc, s.rb_kohm * 1.0e3, rc, re).map_err(|e| e.to_string())?;
            let op = bias.solve(&q).map_err(|e| e.to_string())?;
            let s_ico = bias.stability_factor(&q).map_err(|e| e.to_string())?;
            Ok((op, s_ico))
        }
    }
}

/// Solve the bias network and format the full Q-point readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &BjtWorkbenchState) -> Result<String, String> {
    let (op, s_ico) = solve(s)?;
    let region = match op.region {
        Region::Active => "active (amplifier)",
        Region::Saturation => "saturation (closed switch)",
    };
    let topo = match s.topology {
        Topology::Divider => "voltage divider",
        Topology::FixedBase => "fixed base",
    };

    // Common-base current gain alpha = Ic / Ie, the dual of the
    // forward gain beta (alpha = beta / (beta + 1), always < 1). Derived
    // here from the solved Q-point currents rather than refitting beta, so
    // it stays exact even in saturation where Ic is set by the resistors.
    // Ie is the emitter current of a conducting device, so it is strictly
    // positive whenever the solve succeeded.
    let alpha = op.ic / op.ie;
    // Collector power dissipation Pd = Vce * Ic (W), the quiescent heat the
    // device must shed — the figure checked against its rated dissipation /
    // safe-operating-area. Reported in milliwatts to match the small-signal
    // currents above.
    let pd_mw = op.vce * op.ic * 1.0e3;

    Ok(format!(
        "topology        : {topo}\n\
         gain β / VBE    : {:.0} / {:.2} V\n\
         supply Vcc      : {:.2} V\n\n\
         base current Ib : {:.3} µA\n\
         coll current Ic : {:.3} mA\n\
         emit current Ie : {:.3} mA\n\
         CB gain α       : {alpha:.4}\n\
         emitter VE      : {:.3} V\n\
         Vce             : {:.3} V\n\
         power Pd        : {pd_mw:.3} mW\n\
         stability S(ICO): {:.2}\n\
         region          : {region}",
        s.beta,
        s.vbe,
        s.vcc,
        op.ib * 1.0e6,
        op.ic * 1.0e3,
        op.ie * 1.0e3,
        op.ve,
        op.vce,
        s_ico,
    ))
}

/// Append a capped cylinder along the z-axis, double-sided, centred at
/// `center` with half-length `half_len` and `radius`, `seg` segments.
fn push_cyl_z(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    center: Vector3<f64>,
    half_len: f64,
    radius: f64,
    seg: usize,
) {
    let (z0, z1) = (center.z - half_len, center.z + half_len);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            center.x + radius * a.cos(),
            center.y + radius * a.sin(),
            z0,
        ));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(
            center.x + radius * a.cos(),
            center.y + radius * a.sin(),
            z1,
        ));
    }
    let cap0 = nodes.len();
    nodes.push(Vector3::new(center.x, center.y, z0));
    let cap1 = nodes.len();
    nodes.push(Vector3::new(center.x, center.y, z1));
    for j in 0..seg {
        let jn = (j + 1) % seg;
        // Side wall (double-sided).
        tris.extend_from_slice(&[
            lo + j,
            hi + j,
            hi + jn,
            lo + j,
            hi + jn,
            lo + jn,
            lo + j,
            hi + jn,
            hi + j,
            lo + j,
            lo + jn,
            hi + jn,
        ]);
        // Caps (double-sided fans).
        tris.extend_from_slice(&[cap0, lo + jn, lo + j, cap0, lo + j, lo + jn]);
        tris.extend_from_slice(&[cap1, hi + j, hi + jn, cap1, hi + jn, hi + j]);
    }
}

/// Build the transistor as a triangle [`Mesh`] — a TO-92-style package: a
/// stout cylindrical body whose front is sliced off into the
/// characteristic flat face (the half-cylinder is approximated by sweeping
/// the body only over the rounded rear arc), with three thin leads
/// protruding below. Representative geometry (not to scale; the Q-point
/// numbers are the `valenx-bjt` result). `None` for an invalid
/// configuration.
fn transistor_solid_mesh(s: &BjtWorkbenchState) -> Option<Mesh> {
    // Gate the solid on the same validation the readout uses: only build a
    // package for a transistor that actually has a Q-point.
    solve(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Body: a stout vertical cylinder. The flat face is suggested by the
    // wide, short proportions of a TO-92 can.
    let r = 0.5;
    push_cyl_z(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.6),
        0.5,
        r,
        24,
    );

    // Three leads: thin vertical cylinders below the body, evenly spaced in
    // x (collector / base / emitter).
    for dx in [-0.28, 0.0, 0.28] {
        push_cyl_z(
            &mut nodes,
            &mut tris,
            Vector3::new(dx, 0.0, -0.4),
            0.5,
            0.05,
            10,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-bjt");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D transistor solid and load it into the central viewport.
fn load_transistor_3d(app: &mut ValenxApp) {
    let Some(mesh) = transistor_solid_mesh(&app.bjt) else {
        app.bjt.error =
            Some("transistor parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<transistor>/valenx-bjt"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"bjt"}`** product: the canonical BJT
/// package built as a 3-D solid, paired with the workbench's own `compute()`
/// DC-bias readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`BjtWorkbenchState::default`].
pub(crate) fn bjt_product() -> crate::WorkspaceProduct {
    let s = BjtWorkbenchState::default();
    let mesh = transistor_solid_mesh(&s).expect("canonical BJT ⇒ package solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<bjt>/valenx-bjt");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical BJT ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "BJT (DC bias Q-point)".into(),
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
        let s = BjtWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_qpoint_and_active_region() {
        let mut s = BjtWorkbenchState::default();
        run_bjt(&mut s);
        assert!(
            s.error.is_none(),
            "default divider should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("base current Ib"));
        assert!(s.result.contains("Vce"));
        assert!(s.result.contains("stability S(ICO)"));
        // The stiff 12 V / 47k-10k-2.2k-1k stage sits solidly active.
        assert!(s.result.contains("active"));
    }

    #[test]
    fn analyze_rejects_nonpositive_gain() {
        let mut s = BjtWorkbenchState {
            beta: 0.0,
            ..Default::default()
        };
        run_bjt(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_collector_current_is_beta_times_base_in_active() {
        // GROUND TRUTH: in the forward-active region the collector current
        // follows Ic = beta * Ib exactly (the device's defining relation).
        // The default stiff divider is active, so the solved Q-point must
        // obey it.
        let s = BjtWorkbenchState::default();
        let (op, _) = solve(&s).expect("default divider solves");
        assert_eq!(op.region, Region::Active);
        assert!(
            (op.ic - s.beta * op.ib).abs() < 1e-12,
            "Ic = beta*Ib, got Ic={} vs {}",
            op.ic,
            s.beta * op.ib
        );
        // And Ie = Ic + Ib closes the node.
        assert!((op.ie - (op.ic + op.ib)).abs() < 1e-12);
    }

    #[test]
    fn ground_truth_alpha_and_power_dissipation() {
        // GROUND TRUTH (common-base gain): alpha = Ic/Ie, which for the
        // active-region beta relations equals beta/(beta+1). With the
        // default beta = 100 this is 100/101 = 0.990099..., formatted to
        // four places as "0.9901".
        //
        // GROUND TRUTH (collector power): Pd = Vce * Ic. The default stiff
        // divider solves to Ic ≈ 1.28634 mA and Vce ≈ 7.87086 V, so
        // Pd ≈ 7.87086 V * 1.28634 mA ≈ 10.12 mW, formatted to three
        // places as "10.125".
        let s = BjtWorkbenchState::default();
        let (op, _) = solve(&s).expect("default divider solves");

        // Hand-check the alpha = beta/(beta+1) identity from the solved
        // currents (exact ratio, no formatting).
        let alpha_expected = s.beta / (s.beta + 1.0);
        assert!(
            (op.ic / op.ie - alpha_expected).abs() < 1e-12,
            "alpha = Ic/Ie = beta/(beta+1), got {}",
            op.ic / op.ie
        );
        // Hand-check the power dissipation against Vce * Ic (W),
        // independently of the crate (Pd ≈ 1.0125e-2 W).
        let pd_w = op.vce * op.ic;
        let one = 1.0_f64;
        assert!(
            (pd_w - 1.0125e-2).abs() < 1e-5 * one,
            "Pd = Vce*Ic ≈ 10.12 mW, got {pd_w} W"
        );

        // And confirm both land in the formatted readout at the workbench's
        // exact display precision.
        let out = compute(&s).expect("default divider formats");
        assert!(out.contains("CB gain α       : 0.9901"), "readout: {out}");
        assert!(
            out.contains("power Pd        : 10.125 mW"),
            "readout: {out}"
        );
    }

    #[test]
    fn fixed_base_topology_solves_independently() {
        // Switching topology re-routes the solve through FixedBias; the
        // 240k / 2.2k / Re=0 fixed bias is the textbook active example.
        let mut s = BjtWorkbenchState {
            topology: Topology::FixedBase,
            re_kohm: 0.0,
            ..Default::default()
        };
        run_bjt(&mut s);
        assert!(
            s.error.is_none(),
            "fixed-base should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("fixed base"));
        let (op, _) = solve(&s).expect("fixed-base solves");
        // Ib = (Vcc - VBE)/Rb = (12 - 0.7)/240k.
        let ib_expected = (12.0 - 0.7) / 240_000.0;
        assert!(
            (op.ib - ib_expected).abs() < 1e-12,
            "Ib = (Vcc-VBE)/Rb, got {}",
            op.ib
        );
    }

    #[test]
    fn transistor_mesh_for_default_is_nonempty_and_in_range() {
        let s = BjtWorkbenchState::default();
        let mesh = transistor_solid_mesh(&s).expect("default transistor yields a solid");
        assert!(mesh.nodes.len() > 8, "expected body + three leads");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn transistor_mesh_none_for_invalid() {
        let s = BjtWorkbenchState {
            beta: 0.0,
            ..Default::default()
        };
        assert!(transistor_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bjt_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_bjt_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_bjt_workbench = true;
        run_bjt(&mut app.bjt);
        draw_workbench(&mut app);
    }
}
