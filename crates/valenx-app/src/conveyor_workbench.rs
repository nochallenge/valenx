//! The right-side **Conveyor Workbench** panel — native belt-conveyor
//! throughput and drive-power analysis over `valenx-conveyor`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_conveyor_workbench`,
//! toggled from the View menu. The form sets a belt (bulk density, load
//! area, speed) and an inclined drive path (length, incline angle,
//! friction resistance); "Analyze" computes the mass flow `rho * A * v`,
//! the rated capacity (t/h), the lift power `mdot * g * H` and the total
//! drive power `F * v` as a friction-plus-lift breakdown, and "Show 3-D"
//! loads a representative inclined belt slab on two pulleys into the
//! central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_conveyor::{incline_rise, lift_power, Belt, PowerBreakdown};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Bulk material preset selecting a representative bulk density.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Material {
    /// Dry sand, ~1500 kg/m^3.
    Sand,
    /// Run-of-mine coal, ~850 kg/m^3.
    Coal,
    /// Crushed rock / aggregate, ~1600 kg/m^3.
    Gravel,
    /// Cereal grain, ~770 kg/m^3.
    Grain,
}

impl Material {
    /// Representative bulk density for the preset, kg/m^3.
    fn density(self) -> f64 {
        match self {
            Material::Sand => 1500.0,
            Material::Coal => 850.0,
            Material::Gravel => 1600.0,
            Material::Grain => 770.0,
        }
    }

    /// Short label for the preset.
    fn label(self) -> &'static str {
        match self {
            Material::Sand => "sand",
            Material::Coal => "coal",
            Material::Gravel => "gravel",
            Material::Grain => "grain",
        }
    }
}

/// Persistent form + result state for the Conveyor Workbench.
pub struct ConveyorWorkbenchState {
    /// Selected bulk-material preset (sets the bulk density).
    material: Material,
    /// Material cross-sectional load area `A` (m^2).
    load_area_m2: f64,
    /// Belt travel speed `v` (m/s).
    speed_m_per_s: f64,
    /// Belt path length `L` along the incline (m).
    length_m: f64,
    /// Incline angle above horizontal (degrees).
    incline_deg: f64,
    /// Effective horizontal friction resistance `F_friction` (N).
    friction_n: f64,
    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D belt solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for ConveyorWorkbenchState {
    fn default() -> Self {
        // Sand (1500 kg/m^3) on a 0.20 m^2 load section at 2.0 m/s up a
        // 60 m belt inclined 15 deg, with a 4000 N friction resistance:
        // mdot = 1500 * 0.20 * 2.0 = 600 kg/s, capacity 2160 t/h,
        // H = 60 * sin(15 deg) ~ 15.5 m, lift ~ 91.4 kW.
        Self {
            material: Material::Sand,
            load_area_m2: 0.20,
            speed_m_per_s: 2.0,
            length_m: 60.0,
            incline_deg: 15.0,
            friction_n: 4000.0,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Conveyor Workbench right-side panel. A no-op when the
/// `show_conveyor_workbench` toggle is off.
pub fn draw_conveyor_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_conveyor_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_conveyor_workbench",
        "Conveyor",
        |app, ui| {
            ui.label(
                egui::RichText::new(
                    "native belt-conveyor throughput + drive power · valenx-conveyor",
                )
                .weak()
                .small(),
            );
            ui.separator();

            let s = &mut app.conveyor;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Material").strong());
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.material, Material::Sand, "sand");
                        ui.radio_value(&mut s.material, Material::Coal, "coal");
                        ui.radio_value(&mut s.material, Material::Gravel, "gravel");
                        ui.radio_value(&mut s.material, Material::Grain, "grain");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Belt").strong());
                    ui.horizontal(|ui| {
                        ui.label("load area A (m²)");
                        ui.add(egui::DragValue::new(&mut s.load_area_m2).speed(0.01));
                    });
                    ui.horizontal(|ui| {
                        ui.label("speed v (m/s)");
                        ui.add(egui::DragValue::new(&mut s.speed_m_per_s).speed(0.1));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Drive path").strong());
                    ui.horizontal(|ui| {
                        ui.label("length L (m)");
                        ui.add(egui::DragValue::new(&mut s.length_m).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("incline (deg)");
                        ui.add(egui::DragValue::new(&mut s.incline_deg).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("friction F (N)");
                        ui.add(egui::DragValue::new(&mut s.friction_n).speed(50.0));
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_conveyor(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D").strong())
                        .on_hover_text(
                            "Build a representative inclined belt slab on two pulleys (with a carried material layer) as a 3-D solid and load it into the central viewport to orbit",
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
                        ui.label(egui::RichText::new("Throughput & power").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_conveyor_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.conveyor` borrow is
    // released here): build the belt's 3-D solid and load it.
    if app.conveyor.show_3d_request {
        app.conveyor.show_3d_request = false;
        load_belt_3d(app);
    }
}

/// Validate the form, evaluate the conveyor and format the readout.
fn run_conveyor(s: &mut ConveyorWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Build the validated [`Belt`] from the form, the input both the readout
/// and the 3-D gate need. Extracted so it is unit-testable and shared.
fn belt_of(s: &ConveyorWorkbenchState) -> Result<Belt, String> {
    Belt::new(s.material.density(), s.load_area_m2, s.speed_m_per_s).map_err(|e| e.to_string())
}

/// Evaluate the conveyor and format the full readout, mapping any domain
/// error to a display string. Extracted so it is unit-testable.
fn compute(s: &ConveyorWorkbenchState) -> Result<String, String> {
    let belt = belt_of(s)?;
    let mdot = belt.mass_flow();
    let angle_rad = s.incline_deg.to_radians();
    let rise = incline_rise(s.length_m, angle_rad).map_err(|e| e.to_string())?;
    let lift = lift_power(mdot, rise).map_err(|e| e.to_string())?;
    let power = PowerBreakdown::solve(mdot, s.speed_m_per_s, s.length_m, angle_rad, s.friction_n)
        .map_err(|e| e.to_string())?;
    // Split the resolved drive power into its lift and friction parts.
    // The crate resolves `total` and `lift`; the friction part is the
    // remainder `total - lift` (equal to F_friction * v), and the lift
    // share is what fraction of the motor demand raises material.
    let friction_power = power.total - power.lift;
    let lift_share = if power.total > 0.0 {
        100.0 * power.lift / power.total
    } else {
        0.0
    };

    Ok(format!(
        "material        : {} ({:.0} kg/m³)\n\
         load area / v   : {:.3} m² / {:.2} m/s\n\
         length / incline: {:.1} m / {:.1}°\n\
         friction F      : {:.0} N\n\n\
         mass flow ṁ     : {:.2} kg/s\n\
         volumetric Qv   : {:.4} m³/s\n\
         capacity        : {:.1} t/h\n\
         vertical rise H : {:.3} m\n\
         lift power      : {:.2} kW\n\
         friction power  : {:.2} kW\n\
         lift share      : {:.1} %\n\
         drive tension   : {:.1} N\n\
         drive power     : {:.2} kW",
        s.material.label(),
        s.material.density(),
        s.load_area_m2,
        s.speed_m_per_s,
        s.length_m,
        s.incline_deg,
        s.friction_n,
        mdot,
        belt.volumetric_flow(),
        belt.capacity_tph(),
        rise,
        lift / 1000.0,
        friction_power / 1000.0,
        lift_share,
        power.tension,
        power.total / 1000.0,
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

/// Append an axis-aligned cylinder of radius `r` and half-length `hy`
/// centred at `c`, with its axis along `y` (the belt-width direction),
/// approximated by `SEGMENTS` facets.
fn push_pulley(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    r: f64,
    hy: f64,
) {
    const SEGMENTS: usize = 24;
    let base = nodes.len();
    let step = std::f64::consts::TAU / SEGMENTS as f64;
    for i in 0..SEGMENTS {
        let a = step * i as f64;
        let (sa, ca) = a.sin_cos();
        let dx = r * ca;
        let dz = r * sa;
        nodes.push(c + Vector3::new(dx, -hy, dz));
        nodes.push(c + Vector3::new(dx, hy, dz));
    }
    // Side wall.
    for i in 0..SEGMENTS {
        let i0 = base + 2 * i;
        let i1 = base + 2 * ((i + 1) % SEGMENTS);
        tris.extend_from_slice(&[i0, i1, i1 + 1, i0, i1 + 1, i0 + 1]);
    }
    // Two end caps (fans about each rim's centre node).
    let cap0 = nodes.len();
    nodes.push(c + Vector3::new(0.0, -hy, 0.0));
    let cap1 = nodes.len();
    nodes.push(c + Vector3::new(0.0, hy, 0.0));
    for i in 0..SEGMENTS {
        let a0 = base + 2 * i;
        let a1 = base + 2 * ((i + 1) % SEGMENTS);
        tris.extend_from_slice(&[cap0, a1, a0]);
        tris.extend_from_slice(&[cap1, a0 + 1, a1 + 1]);
    }
}

/// Build the conveyor as a triangle [`Mesh`] — an inclined belt slab on
/// two pulleys with a thin carried-material layer on top. Representative
/// geometry (the belt is tilted by the form's incline angle; lengths are
/// schematic, the throughput / power numbers are the `valenx-conveyor`
/// result). `None` for an invalid configuration.
fn belt_solid_mesh(s: &ConveyorWorkbenchState) -> Option<Mesh> {
    // Gate the geometry on the same validity the analysis enforces.
    belt_of(s).ok()?;
    let angle_rad = s.incline_deg.to_radians();
    incline_rise(s.length_m, angle_rad).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    // Schematic span between pulley centres and a fixed half-width.
    let half_span = 1.0_f64;
    let half_width = 0.35_f64;
    let r = 0.12_f64;

    // Inclined belt slab (long in x, thin in z), tilted about y by the
    // incline angle so the discharge (+x) end rises above the feed end.
    let belt_thick = 0.03_f64;
    let along = Vector3::new(angle_rad.cos(), 0.0, angle_rad.sin());
    let up = Vector3::new(-angle_rad.sin(), 0.0, angle_rad.cos());
    let centre = Vector3::new(0.0, 0.0, r);
    // Eight corners of the tilted slab.
    let base = nodes.len();
    for sx in [-1.0_f64, 1.0] {
        for sy in [-1.0_f64, 1.0] {
            for sz in [-1.0_f64, 1.0] {
                let p = centre
                    + along * (sx * half_span)
                    + Vector3::new(0.0, sy * half_width, 0.0)
                    + up * (sz * belt_thick);
                nodes.push(p);
            }
        }
    }
    // Connect the slab corners (indices base+0..base+7) as a box, with the
    // local (sx,sy,sz) ordering above mapping to bit pattern x*4+y*2+z.
    let faces = [
        [0b000, 0b010, 0b110, 0b100],
        [0b001, 0b101, 0b111, 0b011],
        [0b000, 0b100, 0b101, 0b001],
        [0b010, 0b011, 0b111, 0b110],
        [0b000, 0b001, 0b011, 0b010],
        [0b100, 0b110, 0b111, 0b101],
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

    // Thin carried-material layer riding on the upper (+up) face.
    let mat_centre = centre + up * (belt_thick + 0.02);
    let mat_base = nodes.len();
    for sx in [-1.0_f64, 1.0] {
        for sy in [-1.0_f64, 1.0] {
            for sz in [-1.0_f64, 1.0] {
                let p = mat_centre
                    + along * (sx * half_span * 0.92)
                    + Vector3::new(0.0, sy * half_width * 0.7, 0.0)
                    + up * (sz * 0.025);
                nodes.push(p);
            }
        }
    }
    for f in faces {
        tris.extend_from_slice(&[
            mat_base + f[0],
            mat_base + f[1],
            mat_base + f[2],
            mat_base + f[0],
            mat_base + f[2],
            mat_base + f[3],
        ]);
    }

    // Tail pulley (low, -x end) and head pulley (high, +x end).
    let tail = centre - along * half_span;
    let head = centre + along * half_span;
    push_pulley(&mut nodes, &mut tris, tail, r, half_width);
    push_pulley(&mut nodes, &mut tris, head, r, half_width);

    // Support base under the tail pulley.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(tail.x, 0.0, -0.25),
        Vector3::new(0.2, half_width, 0.15),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-conveyor");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D belt solid and load it into the central viewport.
fn load_belt_3d(app: &mut ValenxApp) {
    let Some(mesh) = belt_solid_mesh(&app.conveyor) else {
        app.conveyor.error =
            Some("conveyor parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<conveyor>/valenx-conveyor"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"conveyor"}`** product: the representative
/// inclined belt slab on two pulleys (with a carried-material layer) built from
/// the canonical sand conveyor (0.20 m² load at 2 m/s up a 60 m belt inclined
/// 15°), paired with the throughput + drive-power readout rows, at a fixed 3/4
/// camera. Registered in [`crate::products_registry`]; the per-tool builder the
/// registry dispatches to. Pure — driven off [`ConveyorWorkbenchState::default`].
pub(crate) fn conveyor_product() -> crate::WorkspaceProduct {
    let s = ConveyorWorkbenchState::default();
    let mesh = belt_solid_mesh(&s).expect("canonical conveyor ⇒ inclined-belt solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<conveyor>/valenx-conveyor");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical conveyor ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Belt conveyor (throughput + power)".into(),
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
        let s = ConveyorWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_capacity_and_power() {
        let mut s = ConveyorWorkbenchState::default();
        run_conveyor(&mut s);
        assert!(
            s.error.is_none(),
            "default conveyor should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("mass flow"));
        assert!(s.result.contains("capacity"));
        assert!(s.result.contains("drive power"));
        // sand: mdot = 1500 * 0.20 * 2.0 = 600 kg/s -> 2160.0 t/h.
        assert!(s.result.contains("600.00 kg/s"));
        assert!(s.result.contains("2160.0 t/h"));
    }

    #[test]
    fn analyze_rejects_zero_load_area() {
        let mut s = ConveyorWorkbenchState {
            load_area_m2: 0.0,
            ..Default::default()
        };
        run_conveyor(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn mass_flow_is_rho_a_v_ground_truth() {
        // Ground truth: mdot = rho * A * v, hand-computed, and capacity
        // is that mass flow times 3.6 (kg/s -> t/h).
        let s = ConveyorWorkbenchState::default();
        let belt = belt_of(&s).expect("default belt is valid");
        let expected = 1500.0 * 0.20 * 2.0; // 600 kg/s
        assert!((belt.mass_flow() - expected).abs() < 1e-9);
        assert!((belt.capacity_tph() - expected * 3.6).abs() < 1e-6);
    }

    #[test]
    fn lift_power_is_mdot_g_h_ground_truth() {
        // Ground truth: P_lift = mdot * g * H with H = L * sin(angle),
        // hand-computed independently of the crate's helpers.
        let s = ConveyorWorkbenchState::default();
        let mdot = belt_of(&s).expect("valid").mass_flow();
        let h: f64 = s.length_m * s.incline_deg.to_radians().sin();
        let expected = mdot * valenx_conveyor::G * h;
        let lift = lift_power(mdot, h).expect("valid");
        assert!((lift - expected).abs() < 1e-6);
    }

    #[test]
    fn friction_power_split_is_total_minus_lift_ground_truth() {
        // Ground truth: the friction part of the resolved drive power is
        // total - lift = F_friction * v, which for the defaults is
        // 4000 N * 2.0 m/s = 8000 W = 8.00 kW exactly (no transcendental,
        // so the substring is libm-independent). The lift share is
        // 100 * lift / total = 91.9 % for the same defaults.
        let mut s = ConveyorWorkbenchState::default();
        run_conveyor(&mut s);
        assert!(
            s.error.is_none(),
            "default conveyor should analyze: {:?}",
            s.error
        );

        // Hand-computed split from the form, independent of the readout.
        let mdot = belt_of(&s).expect("valid").mass_flow();
        let h = s.length_m * s.incline_deg.to_radians().sin();
        let lift = mdot * valenx_conveyor::G * h;
        let expected_friction = s.friction_n * s.speed_m_per_s; // F * v
        assert!(
            (expected_friction - 8000.0).abs() < 1e-9,
            "defaults give F*v = 8000 W, got {expected_friction}"
        );
        // The resolved total is friction + lift; the friction part is the
        // remainder, matching F * v.
        let total = expected_friction + lift;
        assert!(((total - lift) - expected_friction).abs() < 1e-6);

        // The readout surfaces both the friction power and the lift share.
        assert!(
            s.result.contains("friction power  : 8.00 kW"),
            "expected friction power 8.00 kW in:\n{}",
            s.result
        );
        assert!(
            s.result.contains("lift share      : 91.9 %"),
            "expected lift share 91.9 % in:\n{}",
            s.result
        );
    }

    #[test]
    fn belt_mesh_for_default_is_nonempty_and_in_range() {
        let s = ConveyorWorkbenchState::default();
        let mesh = belt_solid_mesh(&s).expect("default conveyor yields a solid");
        assert!(mesh.nodes.len() > 16, "expected slab + material + pulleys");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn belt_mesh_none_for_invalid() {
        let s = ConveyorWorkbenchState {
            speed_m_per_s: 0.0,
            ..Default::default()
        };
        assert!(belt_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_conveyor_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_conveyor_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_conveyor_workbench = true;
        run_conveyor(&mut app.conveyor);
        draw_workbench(&mut app);
    }
}
