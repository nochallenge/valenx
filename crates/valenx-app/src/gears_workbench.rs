//! The right-side **Gears Workbench** panel — native involute-gear design
//! over `valenx-gears`.
//!
//! Mirrors the springs / CFD / FEM workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_gears_workbench`,
//! toggled from the View menu. The form drives a [`valenx_gears::GearSpec`];
//! the "Analyze" button computes the design scalars — circular pitch and the
//! pitch / base / addendum / dedendum diameters — plus the meshing gear ratio
//! against a mating tooth count, and renders them as a monospace readout.

use eframe::egui;

use nalgebra::Vector3;

use valenx_gears::{circular_pitch_mm, full_profile, gear_ratio, GearKind, GearSpec};
use valenx_viz::{project_point, OrbitCamera, ViewDirection};

use crate::ValenxApp;

/// Persistent form + result state for the Gears Workbench.
pub struct GearsWorkbenchState {
    /// Gear family (spur / helical / bevel / worm).
    kind: GearKind,
    /// Module `m` (mm) — pitch diameter = `m × teeth`.
    module_mm: f64,
    /// Tooth count `z`.
    teeth: u32,
    /// Pressure angle (degrees) — standard 20°.
    pressure_angle_deg: f64,
    /// Helix angle (degrees) — 0 for spur, ~20–30 for helical.
    helix_angle_deg: f64,
    /// Face width (mm).
    face_width_mm: f64,
    /// Mating gear's tooth count, for the meshing ratio.
    mate_teeth: u32,
    /// Formatted design readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for GearsWorkbenchState {
    fn default() -> Self {
        // Mirror valenx-gears' standard 20-tooth, 1-module, 20° spur gear.
        Self {
            kind: GearKind::Spur,
            module_mm: 1.0,
            teeth: 20,
            pressure_angle_deg: 20.0,
            helix_angle_deg: 0.0,
            face_width_mm: 10.0,
            mate_teeth: 40,
            result: String::new(),
            error: None,
        }
    }
}

impl GearsWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`
    /// so an agent can discover the name space; each string matches exactly the
    /// caption the form draws (and what the widget is `labelled_by`).
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "module m (mm)",
            "teeth z",
            "pressure angle (°)",
            "helix angle (°)",
            "face width (mm)",
            "mating teeth",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` (the bridge turns it into
    /// a `warn` feed note) — never a panic, and no field is written on error.
    /// Ranges mirror the workbench's own validity checks (module / face width
    /// `> 0`, tooth counts `>= 1`, pressure angle in `(0, 90)`, helix angle in
    /// `[0, 90)`), so a value that would make the design invalid is rejected
    /// here rather than producing a downstream error.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "module m (mm)" => {
                let v = value.as_f64()?;
                if !(v.is_finite() && v > 0.0) {
                    return Err(format!("module m must be > 0, got {v}"));
                }
                self.module_mm = v;
            }
            "teeth z" => {
                let z = value.as_i64()?;
                if !(1..=u32::MAX as i64).contains(&z) {
                    return Err(format!("teeth z must be in [1, {}], got {z}", u32::MAX));
                }
                self.teeth = z as u32;
            }
            "pressure angle (°)" => {
                let v = value.as_f64()?;
                if !(v.is_finite() && v > 0.0 && v < 90.0) {
                    return Err(format!("pressure angle must be in (0, 90), got {v}"));
                }
                self.pressure_angle_deg = v;
            }
            "helix angle (°)" => {
                let v = value.as_f64()?;
                if !(v.is_finite() && (0.0..90.0).contains(&v)) {
                    return Err(format!("helix angle must be in [0, 90), got {v}"));
                }
                self.helix_angle_deg = v;
            }
            "face width (mm)" => {
                let v = value.as_f64()?;
                if !(v.is_finite() && v > 0.0) {
                    return Err(format!("face width must be > 0, got {v}"));
                }
                self.face_width_mm = v;
            }
            "mating teeth" => {
                let z = value.as_i64()?;
                if !(1..=u32::MAX as i64).contains(&z) {
                    return Err(format!(
                        "mating teeth must be in [1, {}], got {z}",
                        u32::MAX
                    ));
                }
                self.mate_teeth = z as u32;
            }
            other => return Err(format!("unknown Gears control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]): the same `Result` string the panel renders
    /// once the gear design has been computed, else the last `error`, else `None`
    /// when it has not been run yet. Read-only — lets an agent read the answer
    /// back after driving a compute, closing the live-driving loop.
    pub fn agent_readout(&self) -> Option<String> {
        if !self.result.is_empty() {
            Some(self.result.clone())
        } else {
            self.error.clone()
        }
    }
}

/// Draw the Gears Workbench right-side panel. A no-op when the
/// `show_gears_workbench` toggle is off.
pub fn draw_gears_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_gears_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_gears_workbench",
        "Gears",
        |app, ui| {
            ui.label(
                egui::RichText::new("native involute-gear design · valenx-gears")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.gears;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Family").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.kind, GearKind::Spur, "Spur");
                        ui.radio_value(&mut s.kind, GearKind::Helical, "Helical");
                        ui.radio_value(&mut s.kind, GearKind::Bevel, "Bevel");
                        ui.radio_value(&mut s.kind, GearKind::Worm, "Worm");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Geometry").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise); the hover text mirrors the
                    // caption for a mouse user.
                    ui.horizontal(|ui| {
                        let m = ui.label("module m (mm)");
                        ui.add(egui::DragValue::new(&mut s.module_mm).speed(0.05))
                            .labelled_by(m.id)
                            .on_hover_text("Module m (mm) — pitch diameter = m·z");
                    });
                    ui.horizontal(|ui| {
                        let z = ui.label("teeth z");
                        ui.add(egui::DragValue::new(&mut s.teeth).speed(1.0))
                            .labelled_by(z.id)
                            .on_hover_text("Tooth count z");
                    });
                    ui.horizontal(|ui| {
                        let pa = ui.label("pressure angle (°)");
                        ui.add(egui::DragValue::new(&mut s.pressure_angle_deg).speed(0.5))
                            .labelled_by(pa.id)
                            .on_hover_text("Pressure angle α (degrees)");
                    });
                    ui.horizontal(|ui| {
                        let ha = ui.label("helix angle (°)");
                        ui.add(egui::DragValue::new(&mut s.helix_angle_deg).speed(0.5))
                            .labelled_by(ha.id)
                            .on_hover_text("Helix angle β (degrees; 0 for spur)");
                    });
                    ui.horizontal(|ui| {
                        let fw = ui.label("face width (mm)");
                        ui.add(egui::DragValue::new(&mut s.face_width_mm).speed(0.5))
                            .labelled_by(fw.id)
                            .on_hover_text("Face width (mm)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Mesh").strong());
                    ui.horizontal(|ui| {
                        let mt = ui.label("mating teeth");
                        ui.add(egui::DragValue::new(&mut s.mate_teeth).speed(1.0))
                            .labelled_by(mt.id)
                            .on_hover_text("Mating gear tooth count (for the ratio)");
                    });

                    // Live hint: the pitch (reference) diameter d = m·z.
                    if s.module_mm > 0.0 {
                        let pd = s.module_mm * s.teeth as f64;
                        ui.label(
                            egui::RichText::new(format!("pitch diameter d ≈ {pd:.2} mm  (m·z)"))
                                .weak()
                                .small(),
                        );
                    }

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_gears(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Design").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }

                    // Live 2-D tooth-profile outline (face-on).
                    if let Some(pts) = preview_profile(s) {
                        ui.add_space(6.0);
                        ui.label(egui::RichText::new("Tooth profile preview").strong());
                        draw_profile_preview(ui, &pts);
                    }
                });
        },
    );
    if close {
        app.show_gears_workbench = false;
    }
}

/// Build the [`GearSpec`] from the form and return the full gear outline as
/// a closed XY polyline (z = 0), best-effort `None` for an invalid spec.
/// Drives the live face-on tooth-profile preview.
fn preview_profile(s: &GearsWorkbenchState) -> Option<Vec<Vector3<f64>>> {
    if !(s.module_mm.is_finite() && s.module_mm > 0.0)
        || s.teeth == 0
        || !(s.pressure_angle_deg.is_finite()
            && s.pressure_angle_deg > 0.0
            && s.pressure_angle_deg < 90.0)
        || !(s.helix_angle_deg.is_finite() && s.helix_angle_deg >= 0.0 && s.helix_angle_deg < 90.0)
        || !(s.face_width_mm.is_finite() && s.face_width_mm > 0.0)
    {
        return None;
    }
    let spec = GearSpec {
        kind: s.kind,
        module_mm: s.module_mm,
        teeth: s.teeth,
        pressure_angle_deg: s.pressure_angle_deg,
        helix_angle_deg: s.helix_angle_deg,
        face_width_mm: s.face_width_mm,
    };
    let outline = full_profile(&spec).ok()?;
    if outline.len() < 2 {
        return None;
    }
    // Lift the 2-D (x, y) outline into z = 0 and close the loop.
    let mut pts: Vec<Vector3<f64>> = outline
        .iter()
        .map(|p| Vector3::new(p[0], p[1], 0.0))
        .collect();
    let first = pts[0];
    pts.push(first);
    Some(pts)
}

/// Draw a closed XY polyline as a face-on (front-view, looking down −Z)
/// wireframe in a fixed-height canvas, the camera auto-framed to the
/// geometry's bounds.
fn draw_profile_preview(ui: &mut egui::Ui, pts: &[Vector3<f64>]) {
    let (response, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 200.0),
        egui::Sense::hover(),
    );
    let rect = response.rect;

    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    for p in pts {
        for k in 0..3 {
            let v = p[k] as f32;
            min[k] = min[k].min(v);
            max[k] = max[k].max(v);
        }
    }

    let mut cam = OrbitCamera::default();
    cam.set_view(ViewDirection::Front);
    cam.frame_bounds(min, max);

    let (w, h) = (rect.width(), rect.height());
    let stroke = egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 200, 255));
    for pair in pts.windows(2) {
        let a = project_point(
            &cam,
            w,
            h,
            [pair[0].x as f32, pair[0].y as f32, pair[0].z as f32],
        );
        let b = project_point(
            &cam,
            w,
            h,
            [pair[1].x as f32, pair[1].y as f32, pair[1].z as f32],
        );
        if let (Some(a), Some(b)) = (a, b) {
            painter.line_segment(
                [
                    egui::pos2(rect.min.x + a.x, rect.min.y + a.y),
                    egui::pos2(rect.min.x + b.x, rect.min.y + b.y),
                ],
                stroke,
            );
        }
    }
}

/// Build a [`GearSpec`] from the form, validate it, and format the
/// design-scalar readout. Extracted from the draw closure so it is
/// unit-testable.
/// Run the gear geometry compute (the in-panel **Analyze** action). Factored
/// out so the button and the product self-test ([`crate::self_test`]) share one
/// path.
pub(crate) fn run(app: &mut ValenxApp) {
    run_gears(&mut app.gears);
}

fn run_gears(s: &mut GearsWorkbenchState) {
    s.error = None;

    if !(s.module_mm.is_finite() && s.module_mm > 0.0) {
        s.error = Some("module must be positive".into());
        return;
    }
    if s.teeth == 0 {
        s.error = Some("tooth count must be at least 1".into());
        return;
    }
    // The base circle d·cos(α) needs 0 < α < 90 (cos α > 0).
    if !(s.pressure_angle_deg.is_finite()
        && s.pressure_angle_deg > 0.0
        && s.pressure_angle_deg < 90.0)
    {
        s.error = Some("pressure angle must be between 0° and 90°".into());
        return;
    }
    if !(s.helix_angle_deg.is_finite() && s.helix_angle_deg >= 0.0 && s.helix_angle_deg < 90.0) {
        s.error = Some("helix angle must be between 0° and 90°".into());
        return;
    }
    if !(s.face_width_mm.is_finite() && s.face_width_mm > 0.0) {
        s.error = Some("face width must be positive".into());
        return;
    }
    if s.mate_teeth == 0 {
        s.error = Some("mating tooth count must be at least 1".into());
        return;
    }

    let spec = GearSpec {
        kind: s.kind,
        module_mm: s.module_mm,
        teeth: s.teeth,
        pressure_angle_deg: s.pressure_angle_deg,
        helix_angle_deg: s.helix_angle_deg,
        face_width_mm: s.face_width_mm,
    };

    let p = circular_pitch_mm(s.module_mm);
    // Convention: this gear is the driver, the mating gear the driven, so the
    // ratio is mate ÷ this (> 1 reduces speed / multiplies torque).
    let ratio = gear_ratio(s.mate_teeth, s.teeth);

    s.result = format!(
        "family        : {}\n\
         module m      : {:.3} mm\n\
         teeth z       : {}\n\
         pressure angle: {:.2}\u{00B0}\n\
         helix angle   : {:.2}\u{00B0}\n\
         face width    : {:.3} mm\n\
         mating teeth  : {}\n\n\
         circular pitch: {:.4} mm  (\u{03C0}\u{00B7}m)\n\
         pitch diameter: {:.3} mm  (m\u{00B7}z)\n\
         base diameter : {:.3} mm  (d\u{00B7}cos \u{03B1})\n\
         addendum dia. : {:.3} mm  (d + 2m)\n\
         dedendum dia. : {:.3} mm  (d \u{2212} 2.5m)\n\
         gear ratio    : {:.3}  (mate \u{00F7} this; >1 reduces speed)",
        s.kind.label(),
        s.module_mm,
        s.teeth,
        s.pressure_angle_deg,
        s.helix_angle_deg,
        s.face_width_mm,
        s.mate_teeth,
        p,
        spec.pitch_diameter_mm(),
        spec.base_diameter_mm(),
        spec.addendum_diameter_mm(),
        spec.dedendum_diameter_mm(),
        ratio,
    );
}

/// Canonical 2-stage spur reducer for the Workbench+Agent **3-D workspace
/// tile**: a ~20:1 speed reducer (1500 rpm / 3 kW input) built from four
/// involute spur gears on two layshafts. Stage 1 meshes an 18-tooth pinion
/// with an 80-tooth gear (4.44:1); stage 2 an 18-tooth pinion with an
/// 81-tooth gear (4.5:1); the overall ratio is `(80/18)·(81/18) ≈ 20.0`.
/// Common geometry: module 2.0 mm, 20° pressure angle, 20 mm face width.
///
/// Each gear is emitted by [`valenx_gears::to_solid_spur`] (a
/// [`valenx_cad::Solid::Mesh`] of involute teeth extruded along +z), then
/// shifted to its shaft centre with [`valenx_cad::Solid::translated`] (which
/// moves the mesh nodes for a mesh-backed solid) and tessellated back to a
/// [`valenx_mesh::Mesh`] via [`valenx_cad::solid_to_mesh`] (a no-op pass-through
/// that returns the cached mesh for mesh-backed solids). The four meshes are
/// merged into one `Tri3` surface — node arrays concatenated and triangle
/// indices re-based by the running node offset — and wrapped as a
/// fully-populated [`crate::types::LoadedMesh`] tagged `<gear>/valenx-2stage`,
/// paired with the design-scalar readout rows. The single source of truth for
/// the agent-bridge gear product (see
/// [`crate::agent_commands::AgentCommand::Show3d`] `kind:"gear"`).
///
/// Self-contained + deterministic. Infallible — the canonical specs are valid
/// (positive face width, ≥3 profile vertices), so every `to_solid_spur` /
/// tessellation succeeds.
pub(crate) fn gear_train_loaded_mesh(
) -> (crate::types::LoadedMesh, Vec<String>, Vec<crate::RigidPart>) {
    use valenx_gears::{contact_ratio, GearSpec};

    const MODULE_MM: f64 = 2.0;
    const PA_DEG: f64 = 20.0;
    const FACE_MM: f64 = 20.0;
    // Input drive conditions (for the readout).
    const INPUT_RPM: f64 = 1500.0;
    const INPUT_KW: f64 = 3.0;
    // Tooth counts: stage-1 pinion→gear, stage-2 pinion→gear.
    const Z1P: u32 = 18;
    const Z1G: u32 = 80;
    const Z2P: u32 = 18;
    const Z2G: u32 = 81;

    let spec = |teeth: u32| GearSpec {
        kind: GearKind::Spur,
        module_mm: MODULE_MM,
        teeth,
        pressure_angle_deg: PA_DEG,
        helix_angle_deg: 0.0,
        face_width_mm: FACE_MM,
    };
    let s1p = spec(Z1P);
    let s1g = spec(Z1G);
    let s2p = spec(Z2P);
    let s2g = spec(Z2G);

    // Centre distances C = m·(zₐ+z_b)/2 keep each meshing pair tangent.
    let c1 = MODULE_MM * (Z1P + Z1G) as f64 / 2.0; // input ↔ layshaft
    let c2 = MODULE_MM * (Z2P + Z2G) as f64 / 2.0; // layshaft ↔ output

    // Two-layshaft layout. Stage 1 sits at z = 0; stage 2 is stacked one face
    // width along +z (a real reducer's gears ride side-by-side on the shaft),
    // so the assembly reads as three parallel shafts:
    //   input  pinion : (x=0,  z=0)
    //   layshaft gear+pinion : (x=c1, z=0) and (x=c1, z=zstack)
    //   output gear : (x=c1+c2, z=zstack)
    let zstack = FACE_MM + 5.0;
    let placements: [(GearSpec, f64, f64); 4] = [
        (s1p.clone(), 0.0, 0.0),
        (s1g.clone(), c1, 0.0),
        (s2p.clone(), c1, zstack),
        (s2g.clone(), c1 + c2, zstack),
    ];

    // Merge the four positioned spur solids into one Tri3 mesh, recording each
    // gear's half-open node range (in placement/concatenation order) so the
    // animation can rotate each gear about its own shaft independently. The
    // pivot is that gear's shaft centre `(cx, 0, cz)`.
    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<u32> = Vec::new();
    let mut gear_ranges: Vec<(std::ops::Range<usize>, [f32; 3])> = Vec::with_capacity(4);
    for (gs, cx, cz) in &placements {
        let solid = valenx_gears::to_solid_spur(gs).expect("canonical spur spec ⇒ solid");
        // Move the mesh-backed solid's nodes to the shaft centre, then take the
        // cached mesh straight back out (solid_to_mesh is a pass-through here).
        let placed = solid.translated(*cx, 0.0, *cz).expect("finite translation");
        let m = valenx_cad::solid_to_mesh(&placed, valenx_cad::DEFAULT_TESS_TOLERANCE)
            .expect("mesh-backed solid ⇒ cached mesh");
        let offset = nodes.len() as u32;
        let start = nodes.len();
        nodes.extend_from_slice(&m.nodes);
        let end = nodes.len();
        gear_ranges.push((start..end, [*cx as f32, 0.0, *cz as f32]));
        for blk in &m.element_blocks {
            if blk.element_type == valenx_mesh::ElementType::Tri3 {
                tris.extend(blk.connectivity.iter().map(|&i| i + offset));
            }
        }
    }

    let mut block = valenx_mesh::ElementBlock::new(valenx_mesh::ElementType::Tri3);
    block.connectivity = tris;
    let mut mesh = valenx_mesh::Mesh::new("valenx-2stage-reducer");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();

    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    let loaded = crate::types::LoadedMesh {
        path: std::path::PathBuf::from("<gear>/valenx-2stage"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    };

    // Design scalars from the real spec API.
    let r1 = gear_ratio(Z1G, Z1P);
    let r2 = gear_ratio(Z2G, Z2P);
    let total = r1 * r2;
    let cr1 = contact_ratio(&s1p, &s1g);
    let cr2 = contact_ratio(&s2p, &s2g);
    let out_rpm = INPUT_RPM / total;
    let lines = vec![
        format!("overall ratio: {total:.2}:1  ({INPUT_RPM:.0} rpm in → {out_rpm:.1} rpm out)"),
        format!("input: {INPUT_RPM:.0} rpm, {INPUT_KW:.0} kW"),
        format!("stage 1: {Z1P}T → {Z1G}T  = {r1:.3}:1"),
        format!("stage 2: {Z2P}T → {Z2G}T  = {r2:.3}:1"),
        format!("module {MODULE_MM:.1} mm · {PA_DEG:.0}° PA · {FACE_MM:.0} mm face"),
        format!("contact ratio: stage 1 {cr1:.2}, stage 2 {cr2:.2}"),
        format!("circular pitch: {:.2} mm", circular_pitch_mm(MODULE_MM)),
    ];

    // Per-gear angular rates for the animation. The TRUE kinematics:
    //   input pinion ω0 ; layshaft gear & pinion share the layshaft = -ω0/r1
    //   (meshing pair counter-rotates ⇒ sign flip); output = -(layshaft)/r2 =
    //   +ω0/(r1·r2). We keep those exact ratios and counter-rotation signs but
    //   PRESENTATION-SCALE the input so the fastest gear reads as ~1.3 rev/s
    //   (≈8 rad/s) instead of the real ~157 rad/s blur. `gear_ranges` is in
    //   placement order: [input pinion, layshaft gear, layshaft pinion, output].
    const PRESENT_W0: f32 = 8.0; // rad/s on the input pinion (~1.27 rev/s)
    let r1f = r1 as f32;
    let r2f = r2 as f32;
    let omegas = [
        PRESENT_W0,               // s1p input pinion
        -PRESENT_W0 / r1f,        // s1g layshaft gear (counter to input)
        -PRESENT_W0 / r1f,        // s2p layshaft pinion (rigid on layshaft)
        PRESENT_W0 / (r1f * r2f), // s2g output (counter to layshaft pinion)
    ];
    let parts: Vec<crate::RigidPart> = gear_ranges
        .into_iter()
        .zip(omegas)
        .map(|((node_range, pivot), rad_per_s)| crate::RigidPart {
            node_range,
            axis: [0.0, 0.0, 1.0], // extrusion axis
            pivot,
            rad_per_s,
        })
        .collect();

    (loaded, lines, parts)
}

/// A fixed 3/4-view [`OrbitCamera`] framing the 2-stage gear-train `mesh`
/// (same `frame_bounds` fit + hero angle as [`crate::rocket_workbench::lv1_camera`]),
/// for the Workbench+Agent gear product's per-tile 3-D view.
pub(crate) fn gear_train_camera(mesh: &valenx_mesh::Mesh) -> OrbitCamera {
    let mut camera = OrbitCamera::default();
    if let Some((min, max)) = crate::mesh_loader::mesh_bounding_box(mesh) {
        camera.frame_bounds(min, max);
    }
    camera.azimuth_deg = 35.0;
    camera.elevation_deg = 22.0;
    camera
}

/// The agent-bridge **`show_3d{kind:"gear"}`** product: the lit 2-stage spur
/// reducer mesh + its ratio readout rows at a fixed 3/4 camera, registered in
/// [`crate::products_registry`]. The per-tool builder the registry dispatches
/// to (so the gear product lives here, not in the shared reducer). Pure.
pub(crate) fn gear_product() -> crate::WorkspaceProduct {
    let (mesh, lines, parts) = gear_train_loaded_mesh();
    let camera = gear_train_camera(&mesh.mesh);
    crate::WorkspaceProduct {
        title: "2-stage spur reducer".into(),
        lines,
        mesh: Some(mesh),
        vertex_colors: None,
        camera,
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        // Animated: each gear counter-rotates about its own shaft so the teeth
        // visibly mesh. Starts paused at t = 0; the tile toolbar drives it.
        animation: Some(crate::ProductAnimation {
            playing: false,
            speed: 1.0,
            t: 0.0,
            motion: crate::ProductMotion::RigidParts(parts),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = GearsWorkbenchState::default();
        // A representative numeric set lands in state.
        s.agent_set("module m (mm)", &AgentValue::Float(2.5))
            .unwrap();
        assert_eq!(s.module_mm, 2.5);
        // An integer caption accepts a whole number.
        s.agent_set("teeth z", &AgentValue::Int(31)).unwrap();
        assert_eq!(s.teeth, 31);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into a numeric field) -> Err.
        assert!(s
            .agent_set("module m (mm)", &AgentValue::Str("big".into()))
            .is_err());
        // Out-of-range numeric -> Err, field unchanged.
        assert!(s.agent_set("teeth z", &AgentValue::Int(0)).is_err());
        assert_eq!(s.teeth, 31, "rejected set leaves the field untouched");
    }

    #[test]
    fn preview_profile_for_default_is_a_closed_outline() {
        let s = GearsWorkbenchState::default();
        let pts = preview_profile(&s).expect("default spec yields a profile");
        assert!(pts.len() >= 4);
        // Closed loop: last point == first.
        assert_eq!(pts.first(), pts.last());
        // Non-zero radial extent + a planar XY profile (z = 0).
        let rmax = pts
            .iter()
            .map(|p| (p.x * p.x + p.y * p.y).sqrt())
            .fold(0.0_f64, f64::max);
        assert!(rmax > 0.0);
        assert!(pts.iter().all(|p| p.z == 0.0));
    }

    #[test]
    fn preview_profile_none_for_invalid_spec() {
        let s = GearsWorkbenchState {
            teeth: 0,
            ..Default::default()
        };
        assert!(preview_profile(&s).is_none());
    }

    #[test]
    fn default_state_is_idle() {
        let s = GearsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_spur_gear() {
        let mut s = GearsWorkbenchState::default();
        run_gears(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        // The readout names the core design scalars.
        assert!(s.result.contains("circular pitch"));
        assert!(s.result.contains("pitch diameter"));
        assert!(s.result.contains("base diameter"));
        assert!(s.result.contains("gear ratio"));
        // Standard m=1, z=20: pitch diameter d = m·z = 20.
        assert!(s.result.contains("20.000 mm"));
        // base = d·cos20° = 20 × 0.9396926 ≈ 18.794.
        assert!(s.result.contains("18.794"));
        // ratio = mate ÷ this = 40 ÷ 20 = 2.
        assert!(s.result.contains("2.000"));
    }

    #[test]
    fn analyze_rejects_zero_teeth() {
        let mut s = GearsWorkbenchState {
            teeth: 0,
            ..Default::default()
        };
        run_gears(&mut s);
        assert!(s.error.is_some());
        assert!(s.result.is_empty());
    }

    #[test]
    fn gear_train_parts_tile_the_mesh_node_count() {
        // The four recorded RigidParts must exactly tile [0, node_count): one
        // contiguous range per gear, in placement order, with no gaps or
        // overlaps — so every node belongs to exactly one rotating body.
        let (mesh, _lines, parts) = gear_train_loaded_mesh();
        assert_eq!(parts.len(), 4, "two-stage reducer fuses four gears");
        let node_count = mesh.mesh.nodes.len();
        assert!(node_count > 0);
        // Ranges are contiguous starting at 0 and ending at node_count.
        assert_eq!(parts[0].node_range.start, 0);
        for w in parts.windows(2) {
            assert_eq!(
                w[0].node_range.end, w[1].node_range.start,
                "gear node ranges must be contiguous (no gap/overlap)"
            );
        }
        assert_eq!(
            parts.last().unwrap().node_range.end,
            node_count,
            "the last gear's range must reach the final node"
        );
        // Each gear owns a non-empty slice.
        for p in &parts {
            assert!(
                p.node_range.end > p.node_range.start,
                "non-empty gear range"
            );
        }
    }

    #[test]
    fn gear_train_parts_counter_rotate_with_true_ratios() {
        // Meshing pairs spin in OPPOSITE directions; the layshaft gear and
        // pinion share a shaft (equal ω). All spin about +z. Speeds follow the
        // true ratios r1, r2 (presentation-scaled but ratio-preserving), so the
        // input is the fastest and the output the slowest.
        let (_mesh, _lines, parts) = gear_train_loaded_mesh();
        // Axis is +z for every gear (the extrusion axis).
        for p in &parts {
            assert_eq!(p.axis, [0.0, 0.0, 1.0]);
        }
        let w_in = parts[0].rad_per_s; // input pinion
        let w_lay_g = parts[1].rad_per_s; // layshaft gear (meshes with input)
        let w_lay_p = parts[2].rad_per_s; // layshaft pinion (rigid on layshaft)
        let w_out = parts[3].rad_per_s; // output gear (meshes with layshaft pinion)
        assert!(w_in > 0.0, "input pinion spins +");
        assert!(w_lay_g < 0.0, "layshaft counter-rotates the input");
        assert!(
            (w_lay_g - w_lay_p).abs() < 1e-6,
            "layshaft gear & pinion share the shaft ⇒ same ω"
        );
        assert!(
            w_out > 0.0,
            "output counter-rotates the layshaft pinion (sign flips back +)"
        );
        // True-ratio magnitudes: |ω_in| > |ω_layshaft| > |ω_out| (the train
        // steps down), and the layshaft ratio matches stage-1 r1.
        assert!(w_in.abs() > w_lay_g.abs());
        assert!(w_lay_p.abs() > w_out.abs());
        let r1 = gear_ratio(80, 18) as f32; // Z1G / Z1P
        assert!(
            (w_in.abs() / w_lay_g.abs() - r1).abs() < 1e-3,
            "input:layshaft speed ratio equals stage-1 gear ratio r1"
        );
        // Presentation scale: the fastest gear is a readable ~6–9 rad/s spin,
        // not the real ~157 rad/s blur.
        assert!(
            (6.0..=9.0).contains(&w_in.abs()),
            "fastest gear presentation-scaled into ~1–1.5 rev/s, got {w_in}"
        );
    }

    #[test]
    fn gear_product_is_animated_rigid_parts() {
        // The agent-bridge gear product carries a paused RigidParts animation
        // with the four gear parts.
        let product = gear_product();
        let anim = product.animation.expect("gear product is animated");
        assert!(!anim.playing, "starts paused at t = 0");
        assert_eq!(anim.t, 0.0);
        assert_eq!(anim.speed, 1.0);
        match anim.motion {
            crate::ProductMotion::RigidParts(parts) => assert_eq!(parts.len(), 4),
            crate::ProductMotion::Turntable { .. } => {
                panic!("gears use per-part rigid motion")
            }
        }
    }

    #[test]
    fn product_animation_clock_math() {
        // Mirror the toolbar's clock advance: t += dt * speed only while playing.
        let mut anim = crate::ProductAnimation {
            playing: true,
            speed: 2.0,
            t: 0.0,
            motion: crate::ProductMotion::Turntable {
                axis: [0.0, 0.0, 1.0],
                pivot: [0.0, 0.0, 0.0],
                rad_per_s: 1.0,
            },
        };
        let dt = 0.5_f32;
        anim.t += dt * anim.speed; // 0.5 * 2.0 = 1.0
        assert!((anim.t - 1.0).abs() < 1e-6);
        // Paused ⇒ the clock would not advance (guarded by `if anim.playing`).
        anim.playing = false;
        let before = anim.t;
        if anim.playing {
            anim.t += dt * anim.speed;
        }
        assert_eq!(anim.t, before, "paused clock holds t");
        // Reset zeroes t.
        anim.t = 0.0;
        assert_eq!(anim.t, 0.0);
    }

    #[test]
    fn analyze_rejects_bad_module_angle_and_mate() {
        for bad in [
            GearsWorkbenchState {
                module_mm: 0.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                pressure_angle_deg: 90.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                face_width_mm: 0.0,
                ..Default::default()
            },
            GearsWorkbenchState {
                mate_teeth: 0,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_gears(&mut s);
            assert!(s.error.is_some());
        }
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
            draw_gears_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_gears_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_gears_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gears_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_gears_workbench = true;
        run_gears(&mut app.gears);
        app.gears.error = Some("invalid gear parameters".to_string());
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The geometry + mesh DragValues are SpinButtons; each must be
        // `labelled_by` its caption (egui clears a DragValue's own Name), so an
        // AI / screen reader can find the control by the caption text.
        let mut app = ValenxApp::default();
        app.show_gears_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // module, teeth, pressure angle, helix angle, face width, mating teeth.
        assert!(
            spin_buttons.len() >= 6,
            "expected the gear numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every gear DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["module m (mm)", "teeth z", "mating teeth"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button and the family radio buttons stay named/invokable.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
        assert!(
            nodes
                .iter()
                .any(|(_, n)| n.role() == Role::RadioButton && n.name() == Some("Spur")),
            "the Spur family radio is a named, selectable node"
        );
    }
}
