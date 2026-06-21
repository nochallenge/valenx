//! The right-side **Creep Workbench** panel — native high-temperature
//! creep and stress-rupture analysis over `valenx-creep`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_creep_workbench`,
//! toggled from the View menu. The form sets an operating temperature and
//! applied stress, a Larson-Miller constant and a master-curve LMP, and a
//! Norton-Bailey secondary-creep law (Arrhenius pre-exponential `A0`,
//! activation energy `Q` and stress exponent `n`). "Analyze" evaluates the
//! Larson-Miller rupture life `t_r = 10^(LMP / T - C)` and the Norton
//! steady-state creep rate `epsilon_dot = A(T) * sigma^n` with its
//! accumulated strain and time-to-strain, and "Show 3-D specimen" loads a
//! representative round-bar creep test specimen into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_creep::larson_miller::rupture_time_hours;
use valenx_creep::norton::NortonLaw;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Creep Workbench.
pub struct CreepWorkbenchState {
    /// Operating (absolute) temperature `T` (kelvin).
    temperature_k: f64,
    /// Applied stress `sigma` (MPa).
    stress_mpa: f64,
    /// Larson-Miller constant `C` (dimensionless; `~20` for many steels).
    lm_constant_c: f64,
    /// Master-curve Larson-Miller parameter `LMP` read at the design
    /// stress (dimensionless, in the K·log10(h) grouping convention).
    lmp: f64,
    /// Norton-Bailey Arrhenius pre-exponential factor `A0`.
    a0: f64,
    /// Norton-Bailey creep activation energy `Q` (J/mol).
    activation_energy_j_per_mol: f64,
    /// Norton-Bailey stress exponent `n`.
    stress_exponent_n: f64,
    /// Service time over which to accumulate secondary-creep strain (h).
    service_time_hours: f64,
    /// Strain limit for the time-to-strain query (e.g. 0.01 = 1 %).
    strain_limit: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D specimen solid (serviced after
    /// the panel draws).
    show_3d_request: bool,
}

impl Default for CreepWorkbenchState {
    fn default() -> Self {
        // A 600 deg C (873 K) component at 100 MPa. With C = 20 and a
        // master-curve LMP of 24000 the Larson-Miller life is ~3.1e7 h;
        // the Norton law A0 = 890.6306, Q = 300 kJ/mol, n = 5 gives
        // A(873 K) ~ 1e-15, a steady-state rate of 1e-5 /h, 1 % strain in
        // 1000 h.
        Self {
            temperature_k: 873.0,
            stress_mpa: 100.0,
            lm_constant_c: 20.0,
            lmp: 24000.0,
            a0: 890.6306,
            activation_energy_j_per_mol: 300000.0,
            stress_exponent_n: 5.0,
            service_time_hours: 1000.0,
            strain_limit: 0.01,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Creep Workbench right-side panel. A no-op when the
/// `show_creep_workbench` toggle is off.
pub fn draw_creep_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_creep_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_creep_workbench",
        "Creep",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Larson-Miller rupture + Norton secondary creep · valenx-creep",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.creep;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Operating point").strong());
                    ui.horizontal(|ui| {
                        ui.label("temperature T (K)");
                        ui.add(egui::DragValue::new(&mut s.temperature_k).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("stress σ (MPa)");
                        ui.add(egui::DragValue::new(&mut s.stress_mpa).speed(1.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Larson-Miller rupture").strong());
                    ui.horizontal(|ui| {
                        ui.label("constant C");
                        ui.add(egui::DragValue::new(&mut s.lm_constant_c).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("master-curve LMP");
                        ui.add(egui::DragValue::new(&mut s.lmp).speed(100.0));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Norton secondary creep").strong());
                    ui.horizontal(|ui| {
                        ui.label("pre-exponential A0");
                        ui.add(egui::DragValue::new(&mut s.a0).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("activation Q (J/mol)");
                        ui.add(egui::DragValue::new(&mut s.activation_energy_j_per_mol).speed(1000.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("stress exponent n");
                        ui.add(egui::DragValue::new(&mut s.stress_exponent_n).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("service time (h)");
                        ui.add(egui::DragValue::new(&mut s.service_time_hours).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("strain limit");
                        ui.add(egui::DragValue::new(&mut s.strain_limit).speed(0.001));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_creep(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D specimen").strong())
                        .on_hover_text(
                            "Build a representative round-bar creep test specimen (a slender gauge section between two wider grip ends) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Creep & rupture").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_creep_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.creep` borrow is
    // released here): build the specimen's 3-D solid and load it.
    if app.creep.show_3d_request {
        app.creep.show_3d_request = false;
        load_specimen_3d(app);
    }
}

/// Validate the form, evaluate the creep models and format the readout.
fn run_creep(s: &mut CreepWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Evaluate the Larson-Miller rupture life and the Norton secondary-creep
/// quantities and format the full readout, mapping any domain error to a
/// display string. Extracted so it is unit-testable.
fn compute(s: &CreepWorkbenchState) -> Result<String, String> {
    // Larson-Miller: invert the master-curve LMP for the life at T.
    let life_hours =
        rupture_time_hours(s.lmp, s.temperature_k, s.lm_constant_c).map_err(|e| e.to_string())?;

    // Norton-Bailey secondary creep with an Arrhenius coefficient A(T).
    let law = NortonLaw::with_arrhenius(
        s.a0,
        s.activation_energy_j_per_mol,
        s.temperature_k,
        s.stress_exponent_n,
    )
    .map_err(|e| e.to_string())?;
    let coefficient = law.coefficient();
    let rate = law.rate_at(s.stress_mpa).map_err(|e| e.to_string())?;
    let strain = law
        .accumulated_strain(s.stress_mpa, s.service_time_hours)
        .map_err(|e| e.to_string())?;
    let time_to_limit = law
        .time_to_strain(s.stress_mpa, s.strain_limit)
        .map_err(|e| e.to_string())?;

    Ok(format!(
        "temperature T   : {:.1} K\n\
         stress σ        : {:.1} MPa\n\
         LM constant C   : {:.2}\n\
         master-curve LMP: {:.0}\n\n\
         rupture life    : {:.4e} h\n\n\
         Norton n        : {:.2}\n\
         coefficient A(T): {:.4e}\n\
         creep rate ε̇    : {:.4e} /h\n\
         strain @ {:.0} h : {:.4e}\n\
         t to {:.3} strain: {:.1} h",
        s.temperature_k,
        s.stress_mpa,
        s.lm_constant_c,
        s.lmp,
        life_hours,
        s.stress_exponent_n,
        coefficient,
        rate,
        s.service_time_hours,
        strain,
        s.strain_limit,
        time_to_limit,
    ))
}

/// Append an outward-facing axial (z-aligned) cylinder of radius `radius`
/// spanning `z0..z1`, centred on the z-axis, to the buffers. `segments` is
/// the number of facets around the circumference; the two end caps are
/// triangle fans. Used to build the round-bar specimen's gauge and grips.
fn push_z_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    radius: f64,
    z0: f64,
    z1: f64,
    segments: usize,
) {
    let base = nodes.len();
    // Bottom ring, top ring, then the two cap centres.
    for &z in &[z0, z1] {
        for k in 0..segments {
            let theta = TAU * (k as f64) / (segments as f64);
            nodes.push(Vector3::new(radius * theta.cos(), radius * theta.sin(), z));
        }
    }
    let bottom_centre = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, z0));
    let top_centre = nodes.len();
    nodes.push(Vector3::new(0.0, 0.0, z1));

    for k in 0..segments {
        let kn = (k + 1) % segments;
        let b0 = base + k;
        let b1 = base + kn;
        let t0 = base + segments + k;
        let t1 = base + segments + kn;
        // Side wall (two triangles per facet, outward-facing).
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
        // Bottom cap (fan, facing -z) and top cap (fan, facing +z).
        tris.extend_from_slice(&[bottom_centre, b1, b0]);
        tris.extend_from_slice(&[top_centre, t0, t1]);
    }
}

/// Build the creep test specimen as a triangle [`Mesh`] — a slender
/// vertical round bar (the gauge section) between two wider cylindrical
/// grip ends. Representative geometry (not to scale; the creep and rupture
/// numbers are the `valenx-creep` result). Gated on the inputs producing a
/// valid Norton law, returning `None` for an invalid configuration.
fn specimen_solid_mesh(s: &CreepWorkbenchState) -> Option<Mesh> {
    // Only build a solid for a physically valid configuration.
    NortonLaw::with_arrhenius(
        s.a0,
        s.activation_energy_j_per_mol,
        s.temperature_k,
        s.stress_exponent_n,
    )
    .ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();
    let seg = 24;

    // Lower grip end (wide), slender gauge section, upper grip end (wide).
    push_z_cylinder(&mut nodes, &mut tris, 0.30, 0.0, 0.40, seg);
    push_z_cylinder(&mut nodes, &mut tris, 0.12, 0.40, 1.10, seg);
    push_z_cylinder(&mut nodes, &mut tris, 0.30, 1.10, 1.50, seg);

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-creep");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D specimen solid and load it into the central viewport.
fn load_specimen_3d(app: &mut ValenxApp) {
    let Some(mesh) = specimen_solid_mesh(&app.creep) else {
        app.creep.error = Some("creep parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<specimen>/valenx-creep"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// Agent-bridge product: the canonical creep workbench as a 3-D solid plus its
/// `compute()` readout rows (see [`crate::products_registry`]).
pub(crate) fn creep_product() -> crate::WorkspaceProduct {
    let s = CreepWorkbenchState::default();
    let mesh = specimen_solid_mesh(&s).expect("canonical creep ⇒ specimen solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<creep>/valenx-specimen");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical creep ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Creep (Larson-Miller)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = CreepWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_rupture_and_creep_rate() {
        let mut s = CreepWorkbenchState::default();
        run_creep(&mut s);
        assert!(
            s.error.is_none(),
            "default specimen should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("rupture life"));
        assert!(s.result.contains("creep rate"));
        // Defaults give A(T) ~ 1e-15, a 1e-5 /h secondary rate and a
        // 1 % strain in exactly 1000 h.
        assert!(s.result.contains("1.0000e-15"));
        assert!(s.result.contains("1.0000e-5"));
        assert!(s.result.contains("1000.0 h"));
    }

    #[test]
    fn analyze_rejects_nonpositive_temperature() {
        let mut s = CreepWorkbenchState {
            temperature_k: 0.0,
            ..Default::default()
        };
        run_creep(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn larson_miller_life_matches_defining_formula() {
        // Ground truth: t_r = 10^(LMP / T - C), evaluated by hand.
        let s = CreepWorkbenchState::default();
        let expected = 10f64.powf(s.lmp / s.temperature_k - s.lm_constant_c);
        let got = rupture_time_hours(s.lmp, s.temperature_k, s.lm_constant_c).unwrap();
        assert!(
            (got - expected).abs() / expected < 1e-12,
            "life mismatch: got {got}, expected {expected}"
        );
    }

    #[test]
    fn norton_rate_ratio_for_doubled_stress_is_two_to_the_n() {
        // Ground truth: epsilon_dot = A sigma^n, so doubling the applied
        // stress multiplies the secondary-creep rate by exactly 2^n.
        let s = CreepWorkbenchState::default();
        let law = NortonLaw::with_arrhenius(
            s.a0,
            s.activation_energy_j_per_mol,
            s.temperature_k,
            s.stress_exponent_n,
        )
        .unwrap();
        let base = law.rate_at(s.stress_mpa).unwrap();
        let doubled = law.rate_at(2.0 * s.stress_mpa).unwrap();
        let ratio = doubled / base;
        assert!(
            (ratio - 2f64.powf(s.stress_exponent_n)).abs() < 1e-6,
            "ratio {ratio} should equal 2^{}",
            s.stress_exponent_n
        );
    }

    #[test]
    fn specimen_mesh_for_default_is_nonempty_and_in_range() {
        let s = CreepWorkbenchState::default();
        let mesh = specimen_solid_mesh(&s).expect("default specimen yields a solid");
        assert!(mesh.nodes.len() > 8, "expected grips + gauge cylinders");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn specimen_mesh_none_for_invalid() {
        let s = CreepWorkbenchState {
            temperature_k: 0.0,
            ..Default::default()
        };
        assert!(specimen_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_creep_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_creep_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_creep_workbench = true;
        run_creep(&mut app.creep);
        draw_workbench(&mut app);
    }
}
