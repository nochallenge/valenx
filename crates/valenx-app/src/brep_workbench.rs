//! The right-side **Part B-Rep (CAD)** workbench — a native, in-house
//! boundary-representation solid modeler built on the [`valenx_truck_cad`]
//! kernel, which sits directly on the `truck` Rust CAD library.
//!
//! The user picks two primitives (`box` / `cylinder` / `sphere`) with
//! labelled, AI-settable size controls, chooses a boolean set operation
//! (`union` / `difference` / `intersection`), and presses **Build**. The
//! kernel builds each primitive as a real closed BRep, applies the
//! boolean (`truck-shapeops`), and tessellates the result into a triangle
//! mesh — the workbench then reports the topology + mesh statistics
//! (solid faces/edges/vertices, mesh triangle/vertex counts, volume,
//! bounding box) and shows a lightweight orthographic preview of the
//! tessellated mesh.
//!
//! It mirrors the other real-time workbenches (`topopt_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_brep_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"brep"`. The
//! bridge can set the controls (`agent_set` / `agent_control_names`),
//! read a status line (`agent_readout`), and fire the build via the
//! RunCommand id `brep.build`.
//!
//! Render depth: this V1 shows the tessellated mesh's statistics plus a
//! simple projected-triangle preview. Wiring the result into the full
//! shared 3-D viewport render path is a follow-up.

use eframe::egui;

use valenx_truck_cad::{BoolOp, BrepError, BrepKernel, Primitive, TriMesh};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Which primitive shape a slot builds. A plain, `Copy` selector enum
/// kept separate from [`valenx_truck_cad::Primitive`] (which also carries
/// the sizes) so the combo box and the size controls stay independent.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PrimKind {
    /// Axis-aligned box.
    Box,
    /// Right circular cylinder.
    Cylinder,
    /// Sphere.
    Sphere,
}

impl PrimKind {
    /// Menu label.
    fn label(self) -> &'static str {
        match self {
            PrimKind::Box => "Box",
            PrimKind::Cylinder => "Cylinder",
            PrimKind::Sphere => "Sphere",
        }
    }

    /// Lowercase id / alias parse for the agent bridge.
    fn from_id(s: &str) -> Option<PrimKind> {
        match s.trim().to_ascii_lowercase().as_str() {
            "box" | "cube" | "block" => Some(PrimKind::Box),
            "cylinder" | "cyl" => Some(PrimKind::Cylinder),
            "sphere" | "ball" => Some(PrimKind::Sphere),
            _ => None,
        }
    }
}

/// The result of a build: the boolean solid's BRep topology counts, the
/// tessellated mesh, and its measured volume.
struct BrepResult {
    /// Tessellated triangle mesh of the boolean result.
    mesh: TriMesh,
    /// BRep topology counts `(faces, edges, vertices)` of the result solid.
    topo: (usize, usize, usize),
    /// Signed volume of the result solid (model units³).
    volume: f64,
    /// Whether the result solid is a closed (watertight) 2-manifold.
    closed: bool,
}

/// Persistent state for the Part B-Rep workbench: the two primitive
/// slots (kind + sizes), the boolean op, and the latest build result.
pub struct BrepWorkbenchState {
    /// Primitive A — shape.
    pub a_kind: PrimKind,
    /// Primitive B — shape.
    pub b_kind: PrimKind,
    /// Boolean op combining A and B.
    pub op: BoolOp,

    /// Box X / cylinder & sphere have their own fields; we keep ONE set
    /// of size fields per slot and interpret them by the slot's kind:
    /// for a box, `(sx, sy, sz)` are the extents; for a cylinder,
    /// `sx = radius`, `sz = height`; for a sphere, `sx = radius`.
    /// Slot A sizes.
    pub a_sx: f64,
    /// Slot A size Y (box only).
    pub a_sy: f64,
    /// Slot A size Z (box / cylinder).
    pub a_sz: f64,
    /// Slot B size X (box dx / cylinder radius / sphere radius).
    pub b_sx: f64,
    /// Slot B size Y (box only).
    pub b_sy: f64,
    /// Slot B size Z (box / cylinder).
    pub b_sz: f64,
    /// Translation applied to slot B before the boolean (so the two
    /// primitives actually overlap). X / Y / Z in model units.
    pub b_off_x: f64,
    /// Slot B offset Y.
    pub b_off_y: f64,
    /// Slot B offset Z.
    pub b_off_z: f64,

    /// Latest build result, or `None` before the first build.
    result: Option<BrepResult>,
    /// Last error string (shown in the panel), cleared on a good build.
    error: Option<String>,
}

impl Default for BrepWorkbenchState {
    fn default() -> Self {
        // A sensible, known-good "punched cube" seed: a unit cube minus a
        // cylinder straddling its mid-face (the documented-robust
        // truck-shapeops input). Build is NOT run until the user (or the
        // bridge) presses Build.
        Self {
            a_kind: PrimKind::Box,
            b_kind: PrimKind::Cylinder,
            op: BoolOp::Difference,
            a_sx: 1.0,
            a_sy: 1.0,
            a_sz: 1.0,
            b_sx: 0.25,
            b_sy: 1.0,
            b_sz: 2.0,
            b_off_x: 0.5,
            b_off_y: 0.5,
            b_off_z: -0.5,
            result: None,
            error: None,
        }
    }
}

impl BrepWorkbenchState {
    /// Assemble a slot's `(kind, sizes)` into a kernel [`Primitive`].
    fn primitive(kind: PrimKind, sx: f64, sy: f64, sz: f64) -> Primitive {
        match kind {
            PrimKind::Box => Primitive::Box {
                dx: sx,
                dy: sy,
                dz: sz,
            },
            PrimKind::Cylinder => Primitive::Cylinder {
                radius: sx,
                height: sz,
            },
            PrimKind::Sphere => Primitive::Sphere { radius: sx },
        }
    }

    /// Build A, build B (translated), apply the boolean, tessellate, and
    /// store the result (or an error). Factored out so the in-panel
    /// **Build** button and the `brep.build` bridge id share one path.
    fn build_now(&mut self) {
        match self.try_build() {
            Ok(res) => {
                self.result = Some(res);
                self.error = None;
            }
            Err(e) => {
                self.error = Some(e.to_string());
                // Keep any previous result on screen so a failed retune
                // doesn't blank the preview; the error line explains why.
            }
        }
    }

    /// The fallible build pipeline. Separated so `build_now` can map the
    /// typed [`BrepError`] into the panel's error line.
    fn try_build(&self) -> Result<BrepResult, BrepError> {
        let kernel = BrepKernel::new();
        let a = kernel.primitive(Self::primitive(
            self.a_kind,
            self.a_sx,
            self.a_sy,
            self.a_sz,
        ))?;
        let b_local = kernel.primitive(Self::primitive(
            self.b_kind,
            self.b_sx,
            self.b_sy,
            self.b_sz,
        ))?;
        // Position B so the two solids overlap. Reuse valenx-cad's rigid
        // translate (it round-trips the truck solid through the kernel
        // wrapper), then unwrap back to a truck solid for the boolean.
        let b = valenx_cad::Solid::from_truck(b_local)
            .translated(self.b_off_x, self.b_off_y, self.b_off_z)
            .map_err(|e| BrepError::Kernel(format!("offset B: {e}")))?;
        let b = match b {
            valenx_cad::Solid::Brep(inner) => inner,
            valenx_cad::Solid::Mesh(_) => {
                return Err(BrepError::Kernel("translate yielded a mesh solid".into()))
            }
        };

        let solid = kernel.boolean(self.op, &a, &b)?;
        let mesh = kernel.tessellate(&solid)?;
        let topo = valenx_truck_cad::topology_counts(&solid);
        let volume = valenx_truck_cad::solid_volume(&solid).unwrap_or(f64::NAN);
        let closed = valenx_truck_cad::is_closed_solid(&solid).unwrap_or(false);
        Ok(BrepResult {
            mesh,
            topo,
            volume,
            closed,
        })
    }

    /// The user-visible captions of every control the agent bridge can
    /// set via `SetControl` (returned by `ListControls`). Order follows
    /// the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Primitive A",
            "A size X",
            "A size Y",
            "A size Z",
            "Primitive B",
            "B size X",
            "B size Y",
            "B size Z",
            "B offset X",
            "B offset Y",
            "B offset Z",
            "Boolean op",
        ]
    }

    /// Set one labelled control by its caption, for the agent
    /// `SetControl` bridge. Fail-loud on an unknown caption / wrong type
    /// / out-of-range; no state is written on error and nothing panics.
    /// Size controls read a positive `f64`; `Primitive A` / `Primitive B`
    /// read a shape id (`box`/`cylinder`/`sphere`); `Boolean op` reads an
    /// op id (`union`/`difference`/`intersection`).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        // Shared positive-size setter for the six size + three offset
        // fields. Offsets may be negative, so only sizes are range-checked.
        fn positive(v: f64, what: &str) -> Result<f64, String> {
            if !v.is_finite() || v <= 0.0 {
                return Err(format!("{what} must be finite and > 0, got {v}"));
            }
            if v > 1.0e4 {
                return Err(format!("{what} must be <= 1e4, got {v}"));
            }
            Ok(v)
        }
        fn finite(v: f64, what: &str) -> Result<f64, String> {
            if !v.is_finite() || v.abs() > 1.0e4 {
                return Err(format!("{what} must be finite and |v| <= 1e4, got {v}"));
            }
            Ok(v)
        }
        match name {
            "Primitive A" => {
                let s = value.as_str()?;
                self.a_kind = PrimKind::from_id(s)
                    .ok_or_else(|| format!("Primitive A must be box/cylinder/sphere, got {s:?}"))?;
            }
            "Primitive B" => {
                let s = value.as_str()?;
                self.b_kind = PrimKind::from_id(s)
                    .ok_or_else(|| format!("Primitive B must be box/cylinder/sphere, got {s:?}"))?;
            }
            "Boolean op" => {
                let s = value.as_str()?;
                self.op = BoolOp::from_id(s).ok_or_else(|| {
                    format!("Boolean op must be union/difference/intersection, got {s:?}")
                })?;
            }
            "A size X" => self.a_sx = positive(value.as_f64()?, "A size X")?,
            "A size Y" => self.a_sy = positive(value.as_f64()?, "A size Y")?,
            "A size Z" => self.a_sz = positive(value.as_f64()?, "A size Z")?,
            "B size X" => self.b_sx = positive(value.as_f64()?, "B size X")?,
            "B size Y" => self.b_sy = positive(value.as_f64()?, "B size Y")?,
            "B size Z" => self.b_sz = positive(value.as_f64()?, "B size Z")?,
            "B offset X" => self.b_off_x = finite(value.as_f64()?, "B offset X")?,
            "B offset Y" => self.b_off_y = finite(value.as_f64()?, "B offset Y")?,
            "B offset Z" => self.b_off_z = finite(value.as_f64()?, "B offset Z")?,
            other => return Err(format!("unknown brep control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the
    /// inputs plus the build outcome (solid topology, mesh counts,
    /// volume, watertight flag) once a build exists. `Some` once built,
    /// `None` before the first build (or after a build error with no
    /// prior result).
    pub fn agent_readout(&self) -> Option<String> {
        if let Some(err) = &self.error {
            // Surface the error as the readout when there's no result to
            // report — fail-loud so the agent sees why a build failed.
            if self.result.is_none() {
                return Some(format!("Part B-Rep build failed: {err}"));
            }
        }
        let r = self.result.as_ref()?;
        let (f, e, v) = r.topo;
        let (min, max) = r.mesh.bounds;
        Some(format!(
            "Part B-Rep \u{00B7} {} {} {} \u{00B7} solid {f}F/{e}E/{v}V \u{00B7} mesh {} tris / {} verts \
             \u{00B7} volume {:.4} \u{00B7} {} \u{00B7} bbox [{:.2},{:.2},{:.2}]..[{:.2},{:.2},{:.2}]",
            self.a_kind.label(),
            self.op.id(),
            self.b_kind.label(),
            r.mesh.triangle_count(),
            r.mesh.vertex_count(),
            r.volume,
            if r.closed { "watertight" } else { "open" },
            min[0], min[1], min[2], max[0], max[1], max[2],
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (build)
// ---------------------------------------------------------------------------

/// Run the build (the in-panel **Build** action). Factored out so the
/// button and the `brep.build` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.brep.build_now();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Part B-Rep workbench. A no-op unless toggled on via
/// View → Part B-Rep (CAD).
pub fn draw_brep_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_brep_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_brep_workbench",
        "Part B-Rep (CAD) \u{2014} truck NURBS / boolean solids",
        brep_workbench_body,
    );
    if close {
        app.show_brep_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn brep_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| brep_workbench_body_inner(app, ui));
}

fn brep_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house B-Rep solid modeling \u{00B7} real NURBS primitives + boolean set-ops on the \
             truck CAD kernel [pick two primitives and a boolean (union / difference / \
             intersection); Build creates closed BRep solids, combines them, and tessellates the \
             result to a triangle mesh].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.brep;

    // --- Primitive A --------------------------------------------------------
    ui.label(egui::RichText::new("Primitive A").strong());
    prim_kind_combo(ui, "brep_a_kind", "Primitive A", &mut s.a_kind);
    egui::Grid::new("brep_a_sizes")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            size_row(ui, "A size X", a_size_x_hint(s.a_kind), &mut s.a_sx);
            if s.a_kind == PrimKind::Box {
                size_row(ui, "A size Y", "Box extent along Y.", &mut s.a_sy);
            }
            if matches!(s.a_kind, PrimKind::Box | PrimKind::Cylinder) {
                size_row(ui, "A size Z", a_size_z_hint(s.a_kind), &mut s.a_sz);
            }
        });

    ui.add_space(6.0);

    // --- Boolean op ---------------------------------------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Boolean op");
        let mut chosen: Option<BoolOp> = None;
        egui::ComboBox::from_id_source("brep_bool_op")
            .selected_text(bool_op_label(s.op))
            .show_ui(ui, |ui| {
                for op in [BoolOp::Union, BoolOp::Difference, BoolOp::Intersection] {
                    if ui
                        .selectable_label(s.op == op, bool_op_label(op))
                        .on_hover_text(bool_op_hint(op))
                        .clicked()
                    {
                        chosen = Some(op);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(op) = chosen {
            s.op = op;
        }
    });

    ui.add_space(6.0);

    // --- Primitive B --------------------------------------------------------
    ui.label(egui::RichText::new("Primitive B").strong());
    prim_kind_combo(ui, "brep_b_kind", "Primitive B", &mut s.b_kind);
    egui::Grid::new("brep_b_sizes")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            size_row(ui, "B size X", a_size_x_hint(s.b_kind), &mut s.b_sx);
            if s.b_kind == PrimKind::Box {
                size_row(ui, "B size Y", "Box extent along Y.", &mut s.b_sy);
            }
            if matches!(s.b_kind, PrimKind::Box | PrimKind::Cylinder) {
                size_row(ui, "B size Z", a_size_z_hint(s.b_kind), &mut s.b_sz);
            }
            // Offsets always apply (position B relative to A).
            offset_row(ui, "B offset X", &mut s.b_off_x);
            offset_row(ui, "B offset Y", &mut s.b_off_y);
            offset_row(ui, "B offset Z", &mut s.b_off_z);
        });

    // --- Build --------------------------------------------------------------
    ui.add_space(6.0);
    if ui
        .button("\u{25B6} Build")
        .on_hover_text(
            "Build both primitives as closed BRep solids, apply the boolean op, and tessellate \
             the result to a triangle mesh.",
        )
        .clicked()
    {
        s.build_now();
    }

    ui.add_space(6.0);
    ui.separator();
    draw_result(s, ui);
}

/// Combo box for a primitive-kind slot, labelled (AI-settable) by
/// `caption`.
fn prim_kind_combo(ui: &mut egui::Ui, id: &str, caption: &str, kind: &mut PrimKind) {
    ui.horizontal(|ui| {
        let lbl = ui.label(caption);
        let mut chosen: Option<PrimKind> = None;
        egui::ComboBox::from_id_source(id)
            .selected_text(kind.label())
            .show_ui(ui, |ui| {
                for k in [PrimKind::Box, PrimKind::Cylinder, PrimKind::Sphere] {
                    if ui.selectable_label(*kind == k, k.label()).clicked() {
                        chosen = Some(k);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(k) = chosen {
            *kind = k;
        }
    });
}

/// A labelled positive-size DragValue row inside a grid.
fn size_row(ui: &mut egui::Ui, caption: &str, hint: &str, value: &mut f64) {
    let lbl = ui.label(caption);
    ui.add(
        egui::DragValue::new(value)
            .speed(0.05)
            .range(0.001..=1.0e4)
            .max_decimals(4),
    )
    .labelled_by(lbl.id)
    .on_hover_text(hint);
    ui.end_row();
}

/// A labelled (sign-allowed) offset DragValue row inside a grid.
fn offset_row(ui: &mut egui::Ui, caption: &str, value: &mut f64) {
    let lbl = ui.label(caption);
    ui.add(
        egui::DragValue::new(value)
            .speed(0.05)
            .range(-1.0e4..=1.0e4)
            .max_decimals(4),
    )
    .labelled_by(lbl.id)
    .on_hover_text("Translate primitive B before the boolean so it overlaps primitive A.");
    ui.end_row();
}

/// Hint for the X-size field, which means different things per kind.
fn a_size_x_hint(kind: PrimKind) -> &'static str {
    match kind {
        PrimKind::Box => "Box extent along X.",
        PrimKind::Cylinder => "Cylinder radius.",
        PrimKind::Sphere => "Sphere radius.",
    }
}

/// Hint for the Z-size field (box height vs cylinder height).
fn a_size_z_hint(kind: PrimKind) -> &'static str {
    match kind {
        PrimKind::Box => "Box extent along Z.",
        PrimKind::Cylinder => "Cylinder height along Z.",
        PrimKind::Sphere => "(unused for a sphere)",
    }
}

fn bool_op_label(op: BoolOp) -> &'static str {
    match op {
        BoolOp::Union => "Union (A \u{222A} B)",
        BoolOp::Difference => "Difference (A \u{2212} B)",
        BoolOp::Intersection => "Intersection (A \u{2229} B)",
    }
}

fn bool_op_hint(op: BoolOp) -> &'static str {
    match op {
        BoolOp::Union => "Weld both solids into one.",
        BoolOp::Difference => "Carve B out of A.",
        BoolOp::Intersection => "Keep only the overlapping region.",
    }
}

// ---------------------------------------------------------------------------
// Result render: stats + orthographic mesh preview
// ---------------------------------------------------------------------------

fn draw_result(s: &BrepWorkbenchState, ui: &mut egui::Ui) {
    if let Some(err) = &s.error {
        ui.label(
            egui::RichText::new(format!("Build error: {err}"))
                .color(egui::Color32::from_rgb(220, 110, 90))
                .strong(),
        );
        ui.add_space(4.0);
    }

    let Some(r) = s.result.as_ref() else {
        ui.label(
            egui::RichText::new(
                "No result yet \u{2014} pick two primitives + a boolean op, then press Build.",
            )
            .italics()
            .weak(),
        );
        return;
    };

    let (f, e, v) = r.topo;
    ui.label(
        egui::RichText::new(format!(
            "Solid: {f} faces \u{00B7} {e} edges \u{00B7} {v} vertices \u{00B7} {}",
            if r.closed { "watertight" } else { "open" },
        ))
        .strong()
        .color(egui::Color32::from_rgb(150, 200, 230)),
    );
    ui.label(format!(
        "Mesh: {} triangles \u{00B7} {} vertices \u{00B7} volume {:.4}",
        r.mesh.triangle_count(),
        r.mesh.vertex_count(),
        r.volume,
    ));

    ui.add_space(6.0);
    ui.label(egui::RichText::new("Preview (orthographic X-Y projection)").strong());
    draw_mesh_preview(&r.mesh, ui);
}

/// A lightweight orthographic preview: project the tessellated triangle
/// mesh onto the X-Y plane and stroke each triangle's edges. Not the full
/// shaded 3-D viewport (that's the render-depth follow-up) — just enough
/// to confirm the boolean produced the expected shape.
fn draw_mesh_preview(mesh: &TriMesh, ui: &mut egui::Ui) {
    let available = ui.available_size();
    let w = available.x.clamp(220.0, 680.0);
    let h = 260.0_f32;
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, h), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(16, 18, 24));

    let (min, max) = mesh.bounds;
    let (dx, dy) = ((max[0] - min[0]) as f32, (max[1] - min[1]) as f32);
    if mesh.positions.is_empty() || dx <= 0.0 && dy <= 0.0 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "empty mesh",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Fit the X-Y bounds into the rect with a margin, preserving aspect.
    let inner = rect.shrink(10.0);
    let span = dx.max(dy).max(1e-6);
    let scale = (inner.width().min(inner.height())) / span;
    let project = |p: &[f64; 3]| -> egui::Pos2 {
        let px = inner.left() + ((p[0] - min[0]) as f32) * scale;
        // Flip Y so +Y is up on screen.
        let py = inner.bottom() - ((p[1] - min[1]) as f32) * scale;
        egui::pos2(px, py)
    };

    let stroke = egui::Stroke::new(0.7, egui::Color32::from_rgb(120, 200, 255));
    // Cap the number of triangles drawn so a very dense mesh stays cheap.
    let max_tris = 4000usize;
    let n = mesh.triangle_count().min(max_tris);
    for t in 0..n {
        let i0 = mesh.triangles[t * 3] as usize;
        let i1 = mesh.triangles[t * 3 + 1] as usize;
        let i2 = mesh.triangles[t * 3 + 2] as usize;
        if i0 >= mesh.positions.len() || i1 >= mesh.positions.len() || i2 >= mesh.positions.len() {
            continue;
        }
        let a = project(&mesh.positions[i0]);
        let b = project(&mesh.positions[i1]);
        let c = project(&mesh.positions[i2]);
        painter.line_segment([a, b], stroke);
        painter.line_segment([b, c], stroke);
        painter.line_segment([c, a], stroke);
    }

    painter.rect_stroke(
        rect,
        2.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(50, 60, 80)),
    );
    let shown = if mesh.triangle_count() > max_tris {
        format!("{} of {} tris", n, mesh.triangle_count())
    } else {
        format!("{} tris", mesh.triangle_count())
    };
    painter.text(
        egui::pos2(rect.right() - 4.0, rect.bottom() - 6.0),
        egui::Align2::RIGHT_BOTTOM,
        shown,
        egui::FontId::monospace(9.0),
        egui::Color32::from_gray(170),
    );
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_seed_is_a_punched_cube() {
        let s = BrepWorkbenchState::default();
        assert_eq!(s.a_kind, PrimKind::Box);
        assert_eq!(s.b_kind, PrimKind::Cylinder);
        assert_eq!(s.op, BoolOp::Difference);
        assert!(s.result.is_none());
    }

    #[test]
    fn build_produces_a_result_and_readout() {
        let mut s = BrepWorkbenchState::default();
        s.build_now();
        assert!(
            s.error.is_none(),
            "default seed should build: {:?}",
            s.error
        );
        let r = s.result.as_ref().expect("a result after build");
        // A punched cube: more than the cube's 6 faces, non-empty mesh.
        assert!(r.topo.0 > 6, "punched cube should have >6 faces");
        assert!(r.mesh.triangle_count() > 0, "mesh has triangles");
        assert!(r.volume > 0.0 && r.volume < 1.0, "punched volume in (0,1)");
        let readout = s.agent_readout().expect("readout after build");
        assert!(readout.contains("Part B-Rep"), "readout: {readout}");
        assert!(
            readout.contains("tris"),
            "readout names mesh tris: {readout}"
        );
    }

    #[test]
    fn agent_set_numeric_and_enum_controls() {
        use crate::agent_commands::AgentValue;
        let mut s = BrepWorkbenchState::default();
        s.agent_set("A size X", &AgentValue::Float(2.0))
            .expect("ax");
        assert!((s.a_sx - 2.0).abs() < 1e-9);
        s.agent_set("B offset Z", &AgentValue::Float(-0.25))
            .expect("boz");
        assert!((s.b_off_z + 0.25).abs() < 1e-9);
        s.agent_set("Primitive A", &AgentValue::Str("sphere".into()))
            .expect("prim a");
        assert_eq!(s.a_kind, PrimKind::Sphere);
        s.agent_set("Boolean op", &AgentValue::Str("union".into()))
            .expect("op");
        assert_eq!(s.op, BoolOp::Union);
    }

    #[test]
    fn agent_set_rejects_bad_values() {
        use crate::agent_commands::AgentValue;
        let mut s = BrepWorkbenchState::default();
        assert!(s.agent_set("A size X", &AgentValue::Float(-1.0)).is_err());
        assert!(s
            .agent_set("A size X", &AgentValue::Float(f64::NAN))
            .is_err());
        assert!(s
            .agent_set("Primitive A", &AgentValue::Str("banana".into()))
            .is_err());
        assert!(s
            .agent_set("Boolean op", &AgentValue::Str("nope".into()))
            .is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn readout_is_none_before_build() {
        let s = BrepWorkbenchState::default();
        assert!(s.agent_readout().is_none(), "no readout before a build");
    }

    #[test]
    fn run_bridge_helper_builds_through_app() {
        let mut app = ValenxApp::default();
        run(&mut app);
        assert!(
            app.brep.result.is_some(),
            "the brep.build bridge helper should produce a result"
        );
    }

    #[test]
    fn control_names_are_listed() {
        let names = BrepWorkbenchState::agent_control_names();
        for c in [
            "Primitive A",
            "A size X",
            "Primitive B",
            "B size X",
            "B offset X",
            "Boolean op",
        ] {
            assert!(names.contains(&c), "missing control name {c}");
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_brep_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_brep_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_brep_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_brep_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by`
        // its caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_brep_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Default seed (box A, cylinder B): A has X/Y/Z (3), B has X/Z (2)
        // plus 3 offsets = at least 8 numeric spin buttons.
        assert!(
            spin_buttons.len() >= 8,
            "expected the numeric size/offset controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        for caption in [
            "A size X",
            "A size Y",
            "A size Z",
            "B size X",
            "B size Z",
            "B offset X",
            "B offset Y",
            "B offset Z",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
