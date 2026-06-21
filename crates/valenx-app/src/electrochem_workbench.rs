//! The right-side **Electrochemistry Workbench** panel — native Nernst
//! cell-potential and Faraday electrolysis analysis over `valenx-electrochem`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_electrochem_workbench`,
//! toggled from the View menu. The form sets a two-electrode cell (cathode and
//! anode standard reduction potentials, electrons transferred `n`, temperature,
//! and reaction quotient `Q`) and an electrolysis deposition (current, time,
//! molar mass, electrons); "Analyze" evaluates the Nernst cell potential
//! `E = E0 - (R T / (n F)) ln Q`, the spontaneity and equilibrium constant, and
//! the Faraday mass `m = (Q M) / (n F)`, and "Show 3-D cell" loads a
//! representative electrochemical-cell solid (an electrolyte tank with two
//! electrode plates dipping in) into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_electrochem::cell::{Cell, Spontaneity};
use valenx_electrochem::faraday::Electrolysis;
use valenx_electrochem::nernst::nernst_potential;
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Electrochemistry Workbench.
pub struct ElectrochemWorkbenchState {
    /// Cathode standard reduction potential `E_cathode` (V vs SHE).
    e_cathode_v: f64,
    /// Anode standard reduction potential `E_anode` (V vs SHE).
    e_anode_v: f64,
    /// Electrons transferred in the balanced overall reaction, `n`.
    electrons: f64,
    /// Absolute temperature `T` (K).
    temperature_k: f64,
    /// Reaction quotient `Q` (dimensionless, reduced-over-oxidised).
    quotient: f64,
    /// Electrolysis current `I` (A).
    current_a: f64,
    /// Electrolysis time `t` (s).
    seconds: f64,
    /// Deposited-substance molar mass `M` (g/mol).
    molar_mass_g_per_mol: f64,
    /// Electrons transferred per formula unit in electrolysis, `n`.
    electrolysis_electrons: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D cell solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ElectrochemWorkbenchState {
    fn default() -> Self {
        // Daniell cell: cathode Cu2+/Cu (+0.34 V), anode Zn2+/Zn (-0.76 V),
        // n = 2 electrons at 298.15 K. At Q = 1 the Nernst term vanishes so
        // E = E0 = 1.10 V (K ~ 10^37, strongly spontaneous). The electrolysis
        // side plates copper: 2.00 A for 1.00 h through Cu2+ (M = 63.546 g/mol,
        // n = 2) deposits ~2.371 g.
        Self {
            e_cathode_v: 0.34,
            e_anode_v: -0.76,
            electrons: 2.0,
            temperature_k: 298.15,
            quotient: 1.0,
            current_a: 2.0,
            seconds: 3600.0,
            molar_mass_g_per_mol: 63.546,
            electrolysis_electrons: 2.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Electrochemistry Workbench right-side panel. A no-op when the
/// `show_electrochem_workbench` toggle is off.
pub fn draw_electrochem_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_electrochem_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_electrochem_workbench",
        "Electrochemistry",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native Nernst cell potential + Faraday electrolysis · valenx-electrochem",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.electrochem;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Cell (Nernst)").strong());
                    ui.horizontal(|ui| {
                        ui.label("E cathode (V)");
                        ui.add(egui::DragValue::new(&mut s.e_cathode_v).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("E anode (V)");
                        ui.add(egui::DragValue::new(&mut s.e_anode_v).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("electrons n");
                        ui.add(egui::DragValue::new(&mut s.electrons).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("temperature (K)");
                        ui.add(egui::DragValue::new(&mut s.temperature_k).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("quotient Q");
                        ui.add(egui::DragValue::new(&mut s.quotient).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Electrolysis (Faraday)").strong());
                    ui.horizontal(|ui| {
                        ui.label("current I (A)");
                        ui.add(egui::DragValue::new(&mut s.current_a).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("time t (s)");
                        ui.add(egui::DragValue::new(&mut s.seconds).speed(10.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("molar mass M (g/mol)");
                        ui.add(egui::DragValue::new(&mut s.molar_mass_g_per_mol).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("electrons n");
                        ui.add(egui::DragValue::new(&mut s.electrolysis_electrons).speed(0.1));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_electrochem(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D cell").strong())
                        .on_hover_text(
                            "Build a representative electrochemical cell (an electrolyte tank with two electrode plates dipping in) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Result").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_electrochem_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.electrochem` borrow is
    // released here): build the cell's 3-D solid and load it.
    if app.electrochem.show_3d_request {
        app.electrochem.show_3d_request = false;
        load_cell_3d(app);
    }
}

/// Validate the form, evaluate the cell + electrolysis and format the readout.
fn run_electrochem(s: &mut ElectrochemWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Human-readable label for a [`Spontaneity`] classification.
fn spontaneity_label(sp: Spontaneity) -> &'static str {
    match sp {
        Spontaneity::Spontaneous => "spontaneous (galvanic)",
        Spontaneity::NonSpontaneous => "non-spontaneous (electrolytic)",
        Spontaneity::Equilibrium => "equilibrium",
    }
}

/// Evaluate the cell + electrolysis and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &ElectrochemWorkbenchState) -> Result<String, String> {
    // Nernst side: build the cell, get the standard cell potential, then apply
    // the Nernst equation at the cell level for the operating-point potential.
    let cell = Cell::new(s.e_cathode_v, s.e_anode_v).map_err(|e| e.to_string())?;
    let e_cell_standard = cell.potential();
    let e_cell = nernst_potential(e_cell_standard, s.electrons, s.temperature_k, s.quotient)
        .map_err(|e| e.to_string())?;
    let spont = spontaneity_label(cell.spontaneity());
    let k = cell
        .equilibrium_constant(s.electrons, s.temperature_k)
        .map_err(|e| e.to_string())?;
    let log10_k = k.log10();

    // Faraday side: charge Q = I t, then mass m = (Q M) / (n F) and moles.
    let electrolysis = Electrolysis::new(s.molar_mass_g_per_mol, s.electrolysis_electrons)
        .map_err(|e| e.to_string())?;
    let charge_c = s.current_a * s.seconds;
    let mass_g = electrolysis
        .mass_from_current(s.current_a, s.seconds)
        .map_err(|e| e.to_string())?;
    let moles = electrolysis.moles(charge_c).map_err(|e| e.to_string())?;

    Ok(format!(
        "E cathode / anode: {:.3} / {:.3} V\n\
         electrons n      : {:.2}\n\
         temperature      : {:.2} K\n\
         quotient Q       : {:.4}\n\n\
         E0 cell          : {:.4} V\n\
         E cell (Nernst)  : {:.4} V\n\
         spontaneity      : {spont}\n\
         log10 K          : {:.2}\n\n\
         current / time   : {:.2} A / {:.0} s\n\
         charge Q         : {:.1} C\n\
         molar mass / n   : {:.3} g/mol / {:.2}\n\
         moles deposited  : {:.6} mol\n\
         mass deposited   : {:.4} g",
        s.e_cathode_v,
        s.e_anode_v,
        s.electrons,
        s.temperature_k,
        s.quotient,
        e_cell_standard,
        e_cell,
        log10_k,
        s.current_a,
        s.seconds,
        charge_c,
        s.molar_mass_g_per_mol,
        s.electrolysis_electrons,
        moles,
        mass_g,
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

/// Build the electrochemical cell as a triangle [`Mesh`] — an electrolyte tank
/// (a wide shallow box) with two electrode plates (the cathode and anode)
/// dipping in from above, and a base. Representative geometry (not to scale;
/// the electrochemistry numbers are the `valenx-electrochem` result). `None`
/// for an invalid configuration (the same validation the readout uses).
fn cell_solid_mesh(s: &ElectrochemWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on the same constructors the readout requires, so an
    // invalid cell / electrolysis yields no solid.
    Cell::new(s.e_cathode_v, s.e_anode_v).ok()?;
    Electrolysis::new(s.molar_mass_g_per_mol, s.electrolysis_electrons).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Electrolyte tank (wide, shallow, centred just above the base).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.45),
        Vector3::new(0.8, 0.4, 0.35),
    );
    // Cathode plate (left electrode, dipping in from the top).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(-0.45, 0.0, 0.7),
        Vector3::new(0.04, 0.3, 0.45),
    );
    // Anode plate (right electrode, dipping in from the top).
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.45, 0.0, 0.7),
        Vector3::new(0.04, 0.3, 0.45),
    );
    // Base.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.05),
        Vector3::new(0.9, 0.5, 0.05),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-electrochem");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D cell solid and load it into the central viewport.
fn load_cell_3d(app: &mut ValenxApp) {
    let Some(mesh) = cell_solid_mesh(&app.electrochem) else {
        app.electrochem.error =
            Some("cell parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<cell>/valenx-electrochem"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"electrochem"}`** product: the canonical
/// electrochemical cell built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`ElectrochemWorkbenchState::default`].
pub(crate) fn electrochem_product() -> crate::WorkspaceProduct {
    let s = ElectrochemWorkbenchState::default();
    let mesh = cell_solid_mesh(&s).expect("canonical cell ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<cell>/valenx-electrochem");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical cell ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Electrochemical cell".into(),
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
        let s = ElectrochemWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_cell_potential_and_mass() {
        let mut s = ElectrochemWorkbenchState::default();
        run_electrochem(&mut s);
        assert!(
            s.error.is_none(),
            "default cell should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("E cell (Nernst)"));
        assert!(s.result.contains("mass deposited"));
        assert!(s.result.contains("spontaneity"));
        // Daniell cell at Q=1: E = E0 = 1.1000 V, strongly spontaneous.
        assert!(s.result.contains("1.1000"));
        assert!(s.result.contains("spontaneous (galvanic)"));
        // Copper plating: 2.00 A * 3600 s through Cu2+ deposits ~2.371 g.
        assert!(s.result.contains("2.3710"));
    }

    #[test]
    fn analyze_rejects_zero_electrons() {
        let mut s = ElectrochemWorkbenchState {
            electrons: 0.0,
            ..Default::default()
        };
        run_electrochem(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn potential_equals_e0_cell_at_q_equals_one() {
        // Ground truth: ln(1) = 0, so the Nernst cell potential at Q = 1 is
        // exactly the standard cell potential E0 = E_cathode - E_anode.
        let s = ElectrochemWorkbenchState::default();
        let e_cell_standard = s.e_cathode_v - s.e_anode_v;
        let e_cell = nernst_potential(e_cell_standard, s.electrons, s.temperature_k, 1.0).unwrap();
        assert!(
            (e_cell - e_cell_standard).abs() < 1e-12,
            "at Q=1, E should equal E0={e_cell_standard}, got {e_cell}"
        );
        assert!(
            (e_cell_standard - 1.10).abs() < 1e-12,
            "Daniell E0 is 1.10 V"
        );
    }

    #[test]
    fn faraday_mass_matches_hand_computation() {
        // Ground truth: m = M I t / (n F). For Cu plating (2.00 A, 3600 s,
        // M = 63.546 g/mol, n = 2) this is 63.546*7200/(2*96485.332...) g.
        use valenx_electrochem::constants::FARADAY_C_PER_MOL;
        let s = ElectrochemWorkbenchState::default();
        let electrolysis =
            Electrolysis::new(s.molar_mass_g_per_mol, s.electrolysis_electrons).unwrap();
        let mass = electrolysis
            .mass_from_current(s.current_a, s.seconds)
            .unwrap();
        let expected = s.molar_mass_g_per_mol * s.current_a * s.seconds
            / (s.electrolysis_electrons * FARADAY_C_PER_MOL);
        assert!(
            (mass - expected).abs() < 1e-9,
            "Faraday mass mismatch: {mass} vs {expected}"
        );
        assert!(
            (mass - 2.371).abs() < 1e-3,
            "Cu plating ~2.371 g, got {mass}"
        );
    }

    #[test]
    fn cell_mesh_for_default_is_nonempty_and_in_range() {
        let s = ElectrochemWorkbenchState::default();
        let mesh = cell_solid_mesh(&s).expect("default cell yields a solid");
        assert!(
            mesh.nodes.len() > 8,
            "expected tank + two electrodes + base"
        );
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn cell_mesh_none_for_invalid() {
        let s = ElectrochemWorkbenchState {
            molar_mass_g_per_mol: 0.0,
            ..Default::default()
        };
        assert!(cell_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_electrochem_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_electrochem_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_electrochem_workbench = true;
        run_electrochem(&mut app.electrochem);
        draw_workbench(&mut app);
    }
}
