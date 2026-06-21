//! The right-side **Transformer Workbench** panel — native ideal
//! two-winding transformer relations over `valenx-transformer`.
//!
//! Mirrors the Heat Transfer / Antenna / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_transformer_workbench`,
//! toggled from the View menu. The form specifies a winding (either by the
//! two turn counts or by a direct turns ratio `a = Np/Ns`), a primary
//! voltage and current, an efficiency, and the excitation frequency and
//! peak core flux; "Analyze" reports the ideal secondary voltage and
//! current (`Vs = Vp/a`, `Is = a*Ip`), the apparent power with its
//! efficiency-derated output and loss, and the transformer-EMF-equation
//! voltage on each winding, and "Show 3-D core" loads a representative
//! laminated E-I core with two coil windings into the central viewport.

use std::f64::consts::TAU;
use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_transformer::impedance::reflect_to_primary;
use valenx_transformer::{apparent_power, induced_emf_rms, volts_per_turn, Efficiency, TurnsRatio};

use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// How the ideal turns ratio is specified in the form.
#[derive(Copy, Clone, PartialEq, Eq)]
enum WindingSpec {
    /// Build the ratio from the two winding turn counts,
    /// `a = Np/Ns` (uses `TurnsRatio::from_turns`).
    Turns,
    /// Build the ratio from a directly-entered value `a`
    /// (uses `TurnsRatio::new`).
    Ratio,
}

/// Persistent form + result state for the Transformer Workbench.
pub struct TransformerWorkbenchState {
    /// Whether the ratio comes from the turn counts or a direct value.
    spec: WindingSpec,
    /// Primary winding turn count `Np`.
    turns_primary: f64,
    /// Secondary winding turn count `Ns`.
    turns_secondary: f64,
    /// Direct turns ratio `a = Np/Ns` (used when `spec` is the ratio mode).
    ratio_a: f64,
    /// Primary (RMS) voltage `Vp` (V).
    voltage_primary: f64,
    /// Primary (RMS) current `Ip` (A).
    current_primary: f64,
    /// Excitation frequency `f` (Hz) for the EMF equation.
    frequency_hz: f64,
    /// Peak mutual core flux `Phi_max` (Wb) for the EMF equation.
    peak_flux_wb: f64,
    /// Efficiency `eta` in `(0, 1]` deriving the real output and loss.
    efficiency: f64,
    /// Secondary-side load impedance magnitude `Zs` (Ω), reflected to the
    /// primary as `Zp = a^2 * Zs` (`reflect_to_primary`).
    load_secondary_ohm: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D core solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for TransformerWorkbenchState {
    fn default() -> Self {
        // A 230 -> 23 V step-down distribution-style winding: Np/Ns =
        // 240/24 = a = 10, so Vs = 23 V and Is = 20 A from a 2 A primary.
        // 460 VA in at eta = 0.97 delivers 446.2 W. At 50 Hz a 0.0043 Wb
        // peak core flux develops ~229 V on the 240-turn primary. An 8 Ω
        // secondary load reflects to a^2 * 8 = 800 Ω at the primary.
        Self {
            spec: WindingSpec::Turns,
            turns_primary: 240.0,
            turns_secondary: 24.0,
            ratio_a: 10.0,
            voltage_primary: 230.0,
            current_primary: 2.0,
            frequency_hz: 50.0,
            peak_flux_wb: 0.0043,
            efficiency: 0.97,
            load_secondary_ohm: 8.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Transformer Workbench right-side panel. A no-op when the
/// `show_transformer_workbench` toggle is off.
pub fn draw_transformer_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_transformer_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_transformer_workbench",
        "Transformer",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native ideal two-winding turns-ratio, EMF & efficiency · valenx-transformer",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.transformer;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Winding ratio").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.spec, WindingSpec::Turns, "from turn counts");
                        ui.radio_value(&mut s.spec, WindingSpec::Ratio, "direct ratio a");
                    });
                    ui.horizontal(|ui| {
                        ui.label("primary turns Np");
                        ui.add(egui::DragValue::new(&mut s.turns_primary).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("secondary turns Ns");
                        ui.add(egui::DragValue::new(&mut s.turns_secondary).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("ratio a = Np/Ns");
                        ui.add(egui::DragValue::new(&mut s.ratio_a).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Primary side").strong());
                    ui.horizontal(|ui| {
                        ui.label("voltage Vp (V)");
                        ui.add(egui::DragValue::new(&mut s.voltage_primary).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("current Ip (A)");
                        ui.add(egui::DragValue::new(&mut s.current_primary).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("efficiency η");
                        ui.add(egui::DragValue::new(&mut s.efficiency).speed(0.005));
                    });
                    ui.horizontal(|ui| {
                        ui.label("load Zs (Ω)")
                            .on_hover_text("secondary-side load impedance, reflected to the primary as Zp = a²·Zs");
                        ui.add(egui::DragValue::new(&mut s.load_secondary_ohm).speed(0.5));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("EMF equation").strong());
                    ui.horizontal(|ui| {
                        ui.label("frequency f (Hz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_hz).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("peak flux Φ (Wb)");
                        ui.add(egui::DragValue::new(&mut s.peak_flux_wb).speed(0.0005));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_transformer(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D core").strong())
                        .on_hover_text(
                            "Build a representative laminated E-I core with two coil windings on the centre limb as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Transformer").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_transformer_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.transformer` borrow is
    // released here): build the core's 3-D solid and load it.
    if app.transformer.show_3d_request {
        app.transformer.show_3d_request = false;
        load_core_3d(app);
    }
}

/// Validate the form, evaluate the transformer and format the readout.
fn run_transformer(s: &mut TransformerWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the ideal turns ratio from the form per the selected winding
/// spec — either `TurnsRatio::from_turns` or `TurnsRatio::new`. Shared by
/// `compute` and the 3-D gate so both honour the same validity.
fn turns_ratio(s: &TransformerWorkbenchState) -> Result<TurnsRatio, String> {
    match s.spec {
        WindingSpec::Turns => {
            TurnsRatio::from_turns(s.turns_primary, s.turns_secondary).map_err(|e| e.to_string())
        }
        WindingSpec::Ratio => TurnsRatio::new(s.ratio_a).map_err(|e| e.to_string()),
    }
}

/// Evaluate the transformer and format the full readout, mapping any
/// domain error to a display string. Extracted so it is unit-testable.
fn compute(s: &TransformerWorkbenchState) -> Result<String, String> {
    let ratio = turns_ratio(s)?;
    let a = ratio.ratio();
    let kind = if ratio.is_step_up() {
        "step-up"
    } else if ratio.is_step_down() {
        "step-down"
    } else {
        "isolation"
    };
    let basis = match s.spec {
        WindingSpec::Turns => "from turn counts",
        WindingSpec::Ratio => "direct ratio a",
    };

    // Ideal voltage / current transformation: Vs = Vp/a, Is = a*Ip.
    let vp = s.voltage_primary;
    let ip = s.current_primary;
    let vs = ratio.secondary_voltage(vp).map_err(|e| e.to_string())?;
    let is = ratio.secondary_current(ip).map_err(|e| e.to_string())?;

    // Power conservation and efficiency derating.
    let pin = apparent_power(vp, ip).map_err(|e| e.to_string())?;
    let eta = Efficiency::new(s.efficiency).map_err(|e| e.to_string())?;
    let eta_val = eta.value();
    let pout = eta.output_power(pin).map_err(|e| e.to_string())?;
    let ploss = eta.power_loss(pin).map_err(|e| e.to_string())?;

    // Reflected (referred) load impedance: a secondary load Zs appears at
    // the primary terminals as Zp = a^2 * Zs.
    let zs_load = s.load_secondary_ohm;
    let zp_reflected = reflect_to_primary(&ratio, zs_load).map_err(|e| e.to_string())?;

    // Transformer EMF equation E = sqrt(2)*pi*f*N*Phi on each winding.
    let np = s.turns_primary;
    let ns = s.turns_secondary;
    let f = s.frequency_hz;
    let phi = s.peak_flux_wb;
    let vpt = volts_per_turn(f, phi).map_err(|e| e.to_string())?;
    let emf_primary = induced_emf_rms(f, np, phi).map_err(|e| e.to_string())?;
    let emf_secondary = induced_emf_rms(f, ns, phi).map_err(|e| e.to_string())?;

    Ok(format!(
        "ratio a         : {a:.4}\n\
         ratio basis     : {basis}\n\
         transformer type: {kind}\n\
         primary V / I   : {vp:.2} V / {ip:.2} A\n\
         secondary V     : {vs:.3} V\n\
         secondary I     : {is:.3} A\n\
         input power     : {pin:.2} VA\n\
         efficiency η    : {eta_val:.3}\n\
         output power    : {pout:.2} W\n\
         power loss      : {ploss:.2} W\n\
         load Zs         : {zs_load:.3} Ω\n\
         reflected Zp    : {zp_reflected:.3} Ω  (a²·Zs)\n\n\
         winding turns   : Np {np:.0} / Ns {ns:.0}\n\
         frequency       : {f:.1} Hz\n\
         peak flux Φ     : {phi:.4} Wb\n\
         volts per turn  : {vpt:.4} V\n\
         EMF primary     : {emf_primary:.2} V\n\
         EMF secondary   : {emf_secondary:.2} V"
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

/// Append a closed cylinder (axis along `z`, centred at `c`) of the given
/// `radius`, half-height and segment count — a coil winding. Side wall
/// plus a top and bottom cap fan.
fn push_cylinder(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    radius: f64,
    half_height: f64,
    segments: usize,
) {
    let base = nodes.len();
    // Bottom ring [0..segments) then top ring [segments..2*segments).
    for &z in &[c.z - half_height, c.z + half_height] {
        for k in 0..segments {
            let theta = TAU * (k as f64) / (segments as f64);
            nodes.push(Vector3::new(
                c.x + radius * theta.cos(),
                c.y + radius * theta.sin(),
                z,
            ));
        }
    }
    // Side wall: two triangles per segment.
    for k in 0..segments {
        let k2 = (k + 1) % segments;
        let b0 = base + k;
        let b1 = base + k2;
        let t0 = base + segments + k;
        let t1 = base + segments + k2;
        tris.extend_from_slice(&[b0, b1, t1, b0, t1, t0]);
    }
    // Cap centres and their fans.
    let bottom_c = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, c.z - half_height));
    let top_c = nodes.len();
    nodes.push(Vector3::new(c.x, c.y, c.z + half_height));
    for k in 0..segments {
        let k2 = (k + 1) % segments;
        tris.extend_from_slice(&[bottom_c, base + k2, base + k]);
        tris.extend_from_slice(&[top_c, base + segments + k, base + segments + k2]);
    }
}

/// Build the transformer as a triangle [`Mesh`] — a laminated E-I core
/// (a stack of rectangular window frames: two yokes plus left, right and
/// centre limbs) with two coil-winding cylinders concentric on the centre
/// limb. Representative geometry (not to scale; the electrical numbers are
/// the `valenx-transformer` result). `None` for an invalid configuration.
fn core_solid_mesh(s: &TransformerWorkbenchState) -> Option<Mesh> {
    if compute(s).is_err() {
        return None;
    }

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Laminated core: a stack of thin window frames along the y axis.
    let lam = 6usize;
    let lam_t = 0.04; // half-thickness of one lamination in y
    let pitch = 0.10; // lamination centre-to-centre spacing
    let mid = (lam as f64 - 1.0) * 0.5;
    for i in 0..lam {
        let yc = (i as f64 - mid) * pitch;
        // Left / right / centre vertical limbs.
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(-0.42, yc, 0.5),
            Vector3::new(0.08, lam_t, 0.5),
        );
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.42, yc, 0.5),
            Vector3::new(0.08, lam_t, 0.5),
        );
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, yc, 0.5),
            Vector3::new(0.10, lam_t, 0.5),
        );
        // Top and bottom yokes closing the window.
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, yc, 0.92),
            Vector3::new(0.5, lam_t, 0.08),
        );
        push_box(
            &mut nodes,
            &mut tris,
            Vector3::new(0.0, yc, 0.08),
            Vector3::new(0.5, lam_t, 0.08),
        );
    }

    // Two coil windings concentric on the centre limb (axis along z).
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.5),
        0.20,
        0.34,
        24,
    );
    push_cylinder(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, 0.5),
        0.15,
        0.30,
        24,
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-transformer");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D core solid and load it into the central viewport.
fn load_core_3d(app: &mut ValenxApp) {
    let Some(mesh) = core_solid_mesh(&app.transformer) else {
        app.transformer.error =
            Some("transformer parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<core>/valenx-transformer"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"transformer"}`** product: the canonical
/// two-winding core built as a 3-D solid, paired with the workbench's own
/// `compute()` readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`TransformerWorkbenchState::default`].
pub(crate) fn transformer_product() -> crate::WorkspaceProduct {
    let s = TransformerWorkbenchState::default();
    let mesh = core_solid_mesh(&s).expect("canonical transformer ⇒ core solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<core>/valenx-transformer");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical transformer ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Transformer (core)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = TransformerWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_ratio_voltage_and_power() {
        let mut s = TransformerWorkbenchState::default();
        run_transformer(&mut s);
        assert!(
            s.error.is_none(),
            "default transformer should analyze: {:?}",
            s.error
        );
        // a = 240/24 = 10, step-down.
        assert!(s.result.contains("10.0000"));
        assert!(s.result.contains("step-down"));
        // Vs = 230/10 = 23 V, Is = 10*2 = 20 A.
        assert!(s.result.contains("23.000"));
        assert!(s.result.contains("20.000"));
        // 460 VA in at eta = 0.97 -> 446.2 W out, 13.8 W loss.
        assert!(s.result.contains("446.20"));
        assert!(s.result.contains("13.80"));
        // EMF on the 240-turn primary at 50 Hz, 0.0043 Wb -> ~229.25 V.
        assert!(s.result.contains("229.25"));
        // An 8 Ω secondary load reflects to a²·Zs = 100·8 = 800 Ω.
        assert!(s.result.contains("8.000 Ω"));
        assert!(s.result.contains("800.000 Ω"));
    }

    #[test]
    fn analyze_rejects_zero_secondary_turns() {
        let mut s = TransformerWorkbenchState {
            turns_secondary: 0.0,
            ..Default::default()
        };
        run_transformer(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn ground_truth_voltage_and_emf_match_closed_form() {
        // Ground truth #1: ideal voltage law Vs = Vp * Ns/Np, hand-computed
        // against the crate's secondary_voltage.
        let (np, ns) = (240.0_f64, 24.0_f64);
        let vp = 230.0_f64;
        let a = TurnsRatio::from_turns(np, ns).unwrap();
        let vs_crate = a.secondary_voltage(vp).unwrap();
        let vs_hand = vp * ns / np; // = 23.0 V
        let dv: f64 = vs_crate - vs_hand;
        assert!(dv.abs() < 1e-12, "Vs got {vs_crate}, hand {vs_hand}");

        // Ground truth #2: transformer EMF equation E = 4.4428829... * f *
        // N * Phi, hand-computed against the crate's induced_emf_rms.
        let (f, phi) = (50.0_f64, 0.0043_f64);
        let emf_crate = induced_emf_rms(f, np, phi).unwrap();
        let emf_hand: f64 = 4.442_882_938_158_366_f64 * f * np * phi;
        let de: f64 = emf_crate - emf_hand;
        assert!(de.abs() < 1e-9, "EMF got {emf_crate}, hand {emf_hand}");
    }

    #[test]
    fn ground_truth_reflected_impedance_is_a_squared_times_load() {
        // Reflected (referred) load law: Zp = a^2 * Zs. With a step-up
        // ratio a = 0.5 and a 16 Ω secondary load, hand-computed
        // Zp = 0.5^2 * 16 = 0.25 * 16 = 4.0 Ω — a step-up makes the load
        // look *smaller* from the primary.
        let s = TransformerWorkbenchState {
            spec: WindingSpec::Ratio,
            ratio_a: 0.5,
            load_secondary_ohm: 16.0,
            ..Default::default()
        };
        let out = compute(&s).expect("step-up config analyzes");

        let a = 0.5_f64;
        let zs = 16.0_f64;
        let zp_hand: f64 = a * a * zs; // = 4.0 Ω
        assert!(
            (zp_hand - 4.0_f64).abs() < 1e-12,
            "hand Zp should be 4 Ω, got {zp_hand}"
        );

        // The crate's own reflect_to_primary must agree with the hand value.
        let ratio = TurnsRatio::new(a).unwrap();
        let zp_crate = reflect_to_primary(&ratio, zs).unwrap();
        assert!(
            (zp_crate - zp_hand).abs() < 1e-12,
            "crate Zp {zp_crate} != hand {zp_hand}"
        );

        // And the formatted readout must surface that exact reflected value.
        assert!(
            out.contains("reflected Zp    : 4.000 Ω"),
            "readout missing reflected impedance line:\n{out}"
        );
        // The load echo line should show the 16 Ω secondary load.
        assert!(
            out.contains("load Zs         : 16.000 Ω"),
            "readout missing load line:\n{out}"
        );
    }

    #[test]
    fn core_mesh_for_default_is_nonempty_and_in_range() {
        let s = TransformerWorkbenchState::default();
        let mesh = core_solid_mesh(&s).expect("default core yields a solid");
        assert!(mesh.nodes.len() > 8, "expected laminated core + windings");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn core_mesh_none_for_invalid() {
        // An efficiency above one is rejected, so the 3-D solid is refused.
        let s = TransformerWorkbenchState {
            efficiency: 1.5,
            ..Default::default()
        };
        assert!(core_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_transformer_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_transformer_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_transformer_workbench = true;
        run_transformer(&mut app.transformer);
        draw_workbench(&mut app);
    }
}
