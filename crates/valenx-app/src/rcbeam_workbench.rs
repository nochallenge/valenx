//! The right-side **RC Beam Workbench** panel — native singly-reinforced
//! rectangular reinforced-concrete beam flexure over `valenx-rcbeam`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_rcbeam_workbench`,
//! toggled from the View menu. The form sets the section geometry (width
//! `b`, effective depth `d`), the materials (concrete strength `fc'`, steel
//! yield `fy`) and the tension-steel area `As`; "Analyze" evaluates the
//! Whitney equivalent-stress-block flexure equations and reports the
//! stress-block depth `a`, lever arm `d - a/2`, nominal moment `Mn`, design
//! strength `phi*Mn`, reinforcement ratio `rho`, balanced ratio `rho_b` and
//! the under-reinforced (ductile) check; "Show 3-D beam" loads a
//! representative concrete-beam solid with tension rebars into the central
//! viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;
use valenx_rcbeam::BeamSection;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Steel elastic modulus `Es` (MPa) used for the balanced-ratio /
/// under-reinforced check — the textbook ~200_000 MPa for mild steel.
const STEEL_MODULUS_MPA: f64 = 200_000.0;

/// Stress-block depth factor `beta1 = a/c` used for the balanced-ratio
/// check — ACI-318's `0.85` for `fc' <= 28 MPa.
const BETA1: f64 = 0.85;

/// Persistent form + result state for the RC Beam Workbench.
pub struct RcBeamWorkbenchState {
    /// Section width `b` (mm).
    width_mm: f64,
    /// Effective depth `d` (mm) — extreme compression fibre to tension-steel
    /// centroid.
    effective_depth_mm: f64,
    /// Specified concrete compressive strength `fc'` (MPa).
    fc_prime_mpa: f64,
    /// Specified steel yield strength `fy` (MPa).
    fy_mpa: f64,
    /// Total tension-reinforcement area `As` (mm²).
    area_steel_mm2: f64,
    /// Formatted capacity readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D beam solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for RcBeamWorkbenchState {
    fn default() -> Self {
        // A 300 x 550 mm section, fc' = 30 MPa, fy = 420 MPa, As = 1500 mm^2
        // (SI units: mm, MPa, mm^2 -> N·mm). rho ~ 0.0091 < rho_b ~ 0.0304,
        // so under-reinforced; a ~ 82.35 mm, Mn ~ 320.56 kN·m.
        Self {
            width_mm: 300.0,
            effective_depth_mm: 550.0,
            fc_prime_mpa: 30.0,
            fy_mpa: 420.0,
            area_steel_mm2: 1500.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the RC Beam Workbench right-side panel. A no-op when the
/// `show_rcbeam_workbench` toggle is off.
pub fn draw_rcbeam_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rcbeam_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_rcbeam_workbench",
        "RC Beam",
        |app, ui| {
            ui.label(
                egui::RichText::new("native reinforced-concrete beam flexure · valenx-rcbeam")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.rcbeam;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Section").strong());
                    ui.horizontal(|ui| {
                        ui.label("width b (mm)");
                        ui.add(egui::DragValue::new(&mut s.width_mm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("effective depth d (mm)");
                        ui.add(egui::DragValue::new(&mut s.effective_depth_mm).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Materials").strong());
                    ui.horizontal(|ui| {
                        ui.label("concrete fc' (MPa)");
                        ui.add(egui::DragValue::new(&mut s.fc_prime_mpa).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("steel fy (MPa)");
                        ui.add(egui::DragValue::new(&mut s.fy_mpa).speed(5.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Reinforcement").strong());
                    ui.horizontal(|ui| {
                        ui.label("tension steel As (mm²)");
                        ui.add(egui::DragValue::new(&mut s.area_steel_mm2).speed(10.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_beam(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D beam").strong())
                        .on_hover_text(
                            "Build a representative reinforced-concrete beam (a rectangular concrete prism with tension rebars near the bottom) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Flexural capacity").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_rcbeam_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.rcbeam` borrow is
    // released here): build the beam's 3-D solid and load it.
    if app.rcbeam.show_3d_request {
        app.rcbeam.show_3d_request = false;
        load_beam_3d(app);
    }
}

/// Validate the form, evaluate the section and format the readout.
fn run_beam(s: &mut RcBeamWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the beam section and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &RcBeamWorkbenchState) -> Result<String, String> {
    let section = BeamSection::new(
        s.width_mm,
        s.effective_depth_mm,
        s.fc_prime_mpa,
        s.fy_mpa,
        s.area_steel_mm2,
    )
    .map_err(|e| e.to_string())?;

    let a = section.stress_block_depth();
    let jd = section.lever_arm().map_err(|e| e.to_string())?;
    let mn = section.nominal_moment().map_err(|e| e.to_string())?;
    let phi_mn = section.design_moment_default().map_err(|e| e.to_string())?;
    let rho = section.reinforcement_ratio();
    let rho_b = section
        .balanced_ratio(BETA1, STEEL_MODULUS_MPA)
        .map_err(|e| e.to_string())?;
    let under = section
        .is_under_reinforced(BETA1, STEEL_MODULUS_MPA)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "section b x d   : {:.1} x {:.1} mm\n\
         materials       : fc' {:.1} MPa / fy {:.1} MPa\n\
         tension steel As: {:.1} mm²\n\n\
         stress block a  : {:.2} mm\n\
         lever arm d-a/2 : {:.2} mm\n\
         nominal Mn      : {:.2} kN·m\n\
         design phi*Mn   : {:.2} kN·m  (phi = 0.90)\n\
         steel ratio rho : {:.5}\n\
         balanced rho_b  : {:.5}\n\
         section is      : {}",
        s.width_mm,
        s.effective_depth_mm,
        s.fc_prime_mpa,
        s.fy_mpa,
        s.area_steel_mm2,
        a,
        jd,
        mn / 1.0e6,
        phi_mn / 1.0e6,
        rho,
        rho_b,
        if under {
            "under-reinforced (ductile)"
        } else {
            "NOT under-reinforced"
        },
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

/// Append a (double-sided) cylinder whose axis runs along `+x`, spanning
/// `base.x ..= base.x + length` with circle centre `(base.y, base.z)`.
fn push_cyl_x(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    base: Vector3<f64>,
    length: f64,
    r: f64,
    seg: usize,
) {
    let (x0, x1) = (base.x, base.x + length);
    let lo = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x0, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    let hi = nodes.len();
    for j in 0..seg {
        let a = j as f64 / seg as f64 * TAU;
        nodes.push(Vector3::new(x1, base.y + r * a.cos(), base.z + r * a.sin()));
    }
    for j in 0..seg {
        let jn = (j + 1) % seg;
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
    }
}

/// Build the reinforced-concrete beam as a triangle [`Mesh`] — a long
/// rectangular concrete prism (axis along `+x`) with a few tension rebars
/// running near the bottom face. Representative geometry (not to scale; the
/// flexural numbers are the `valenx-rcbeam` result). `None` for an invalid
/// configuration.
fn beam_solid_mesh(s: &RcBeamWorkbenchState) -> Option<Mesh> {
    // Gate on the real section being constructible *and* its capacity being
    // evaluable (degenerate over-reinforced sections have no valid solid).
    let section = BeamSection::new(
        s.width_mm,
        s.effective_depth_mm,
        s.fc_prime_mpa,
        s.fy_mpa,
        s.area_steel_mm2,
    )
    .ok()?;
    section.nominal_moment().ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Concrete prism: long span (x), section half-width (y) and half-depth
    // (z) scaled from b / d. The beam spans x in [0, span].
    let span = 4.0_f64;
    let half_w = (s.width_mm / 1000.0).max(0.05) * 0.5;
    let half_d = (s.effective_depth_mm / 1000.0).max(0.05) * 0.5;
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(span * 0.5, 0.0, half_d),
        Vector3::new(span * 0.5, half_w, half_d),
    );

    // Three tension rebars running the full span near the bottom face.
    let r_bar = (half_w * 0.12).max(0.01);
    let z_bar = r_bar * 1.6; // sit just above the bottom face
    let n_bars = 3usize;
    for k in 0..n_bars {
        let frac = if n_bars > 1 {
            k as f64 / (n_bars - 1) as f64
        } else {
            0.5
        };
        // Span the bar centres across the width, inset from the faces.
        let y = -half_w * 0.7 + frac * (1.4 * half_w);
        push_cyl_x(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, y, z_bar),
            span,
            r_bar,
            12,
        );
    }

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-rcbeam");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D beam solid and load it into the central viewport.
fn load_beam_3d(app: &mut ValenxApp) {
    let Some(mesh) = beam_solid_mesh(&app.rcbeam) else {
        app.rcbeam.error = Some("beam parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<beam>/valenx-rcbeam"),
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
        let s = RcBeamWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_moment_and_ratios() {
        let mut s = RcBeamWorkbenchState::default();
        run_beam(&mut s);
        assert!(
            s.error.is_none(),
            "default beam should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("nominal Mn"));
        assert!(s.result.contains("design phi*Mn"));
        assert!(s.result.contains("steel ratio rho"));
        assert!(s.result.contains("balanced rho_b"));
        // 300 x 550, fc'=30, fy=420, As=1500 -> Mn ~ 320.56 kN·m, and the
        // section is under-reinforced (rho ~ 0.0091 < rho_b ~ 0.0304).
        assert!(s.result.contains("320.56"), "result: {}", s.result);
        assert!(s.result.contains("under-reinforced (ductile)"));
    }

    #[test]
    fn analyze_rejects_zero_width() {
        let mut s = RcBeamWorkbenchState {
            width_mm: 0.0,
            ..Default::default()
        };
        run_beam(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn capacity_matches_whitney_stress_block_closed_form() {
        // Ground truth (ACI Whitney equivalent stress block, SI units):
        //   a  = As*fy / (0.85*fc'*b)
        //   Mn = As*fy*(d - a/2)
        // For the default section (b=300, d=550, fc'=30, fy=420, As=1500):
        //   a  = 1500*420 / (0.85*30*300) = 82.352941... mm
        //   Mn = 1500*420*(550 - a/2)     = 320_558_823.5... N·mm
        let s = RcBeamWorkbenchState::default();
        let section = BeamSection::new(
            s.width_mm,
            s.effective_depth_mm,
            s.fc_prime_mpa,
            s.fy_mpa,
            s.area_steel_mm2,
        )
        .unwrap();

        let a_expected = 1500.0 * 420.0 / (0.85 * 30.0 * 300.0);
        assert!(
            (section.stress_block_depth() - a_expected).abs() < 1.0e-6,
            "a = {}, expected {a_expected}",
            section.stress_block_depth()
        );

        let mn_expected = 1500.0 * 420.0 * (550.0 - a_expected / 2.0);
        let mn = section.nominal_moment().unwrap();
        assert!(
            (mn - mn_expected).abs() < 1.0e-3,
            "Mn = {mn}, expected {mn_expected}"
        );
        // ~320.56 kN·m.
        assert!(
            (mn / 1.0e6 - 320.558_823_5).abs() < 1.0e-3,
            "Mn = {mn} N·mm"
        );

        // phi*Mn applies the default tension-controlled phi = 0.90.
        let phi_mn = section.design_moment_default().unwrap();
        assert!((phi_mn - 0.90 * mn).abs() < 1.0e-3, "phi*Mn = {phi_mn}");
    }

    #[test]
    fn beam_mesh_for_default_is_nonempty_and_in_range() {
        let s = RcBeamWorkbenchState::default();
        let mesh = beam_solid_mesh(&s).expect("default beam yields a solid");
        assert!(mesh.nodes.len() > 8, "expected concrete prism + rebars");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn beam_mesh_none_for_invalid() {
        let s = RcBeamWorkbenchState {
            width_mm: 0.0,
            ..Default::default()
        };
        assert!(beam_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rcbeam_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_rcbeam_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rcbeam_workbench = true;
        run_beam(&mut app.rcbeam);
        draw_workbench(&mut app);
    }
}
