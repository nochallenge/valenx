//! Panel 8 — **Macromolecular Structure** (`valenx-biostruct`).
//!
//! Load a macromolecular structure (PDB or mmCIF), compute a per-chain
//! secondary-structure / composition summary, classify the
//! Ramachandran φ/ψ distribution, and Kabsch-superpose two structures'
//! Cα atoms for an RMSD — all native `valenx-biostruct` calls.

use eframe::egui;
use nalgebra::Point3;

use valenx_biostruct::analyze::StructureReport;
use valenx_biostruct::geometry::ramachandran::summarize as rama_summarize;
use valenx_biostruct::io::read_structure;
use valenx_biostruct::structure::Structure;
use valenx_biostruct::superpose::{kabsch, rmsd};

use super::common;
use super::molecule_view::{self, ViewMolecule};
use crate::molviz::{self, AtomAttr, BackbonePoint, ColorScheme, MolvizParams, Representation};
use crate::ValenxApp;

/// Which structure sub-tool is showing.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Tool {
    #[default]
    Analyze,
    Ramachandran,
    Superpose,
}

/// Snapshot of every editable input the Biostruct panel owns.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct BiostructSnapshot {
    pub(crate) structure_a: String,
    pub(crate) structure_b: String,
    pub(crate) clash_tolerance: f64,
}

/// Form + result state for the Macromolecular Structure panel.
pub struct BiostructPanel {
    tool: Tool,
    /// Structure text for the analyze / Ramachandran tools, and the
    /// mobile structure of the superpose tool.
    structure_a: String,
    /// Reference structure text for the superpose tool.
    structure_b: String,
    /// Steric-clash tolerance (Å) for the analysis report.
    clash_tolerance: f64,
    error: Option<String>,
    result: String,
    /// Undo / redo over both structure-text inputs + the clash tol.
    history: crate::undo::History<BiostructSnapshot>,
    /// Which molecular-viewer representation the "Show in 3D viewport" button
    /// builds. Reactive picker; default [`Representation::BallAndStick`].
    pub(crate) representation: Representation,
    /// How the viewer colours atoms. Reactive picker; default
    /// [`ColorScheme::Element`] (the CPK palette). When non-default
    /// (chain / residue / B-factor), "Show in 3D viewport" builds a
    /// per-vertex-coloured mesh and uploads the colours so the viewport
    /// renders the scheme instead of monochrome.
    pub(crate) color_scheme: ColorScheme,
    /// Per-representation mesh-generation tunables (surface grid resolution,
    /// cartoon tube, ball/stick radii). Mutated by the picker's sliders.
    pub(crate) molviz_params: MolvizParams,
}

impl BiostructPanel {
    fn snapshot(&self) -> BiostructSnapshot {
        BiostructSnapshot {
            structure_a: self.structure_a.clone(),
            structure_b: self.structure_b.clone(),
            clash_tolerance: self.clash_tolerance,
        }
    }
    fn restore(&mut self, s: BiostructSnapshot) {
        self.structure_a = s.structure_a;
        self.structure_b = s.structure_b;
        self.clash_tolerance = s.clash_tolerance;
    }
    pub fn undo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(prev) = self.history.undo(current) {
            self.restore(prev);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn redo_edit(&mut self) -> bool {
        let current = self.snapshot();
        if let Some(next) = self.history.redo(current) {
            self.restore(next);
            self.error = None;
            true
        } else {
            false
        }
    }
    pub fn can_undo(&self) -> bool {
        self.history.can_undo()
    }
    pub fn can_redo(&self) -> bool {
        self.history.can_redo()
    }
}

/// A minimal 3-residue glycine peptide PDB — enough for the analysis +
/// Ramachandran tools to produce real output without file I/O.
const DEMO_PDB: &str = "\
ATOM      1  N   GLY A   1      -1.204   1.045   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       1.250   0.881   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       1.150   2.100   0.000  1.00  0.00           O
ATOM      5  N   GLY A   2       2.430   0.300   0.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       3.720   0.960   0.000  1.00  0.00           C
ATOM      7  C   GLY A   2       4.880   0.000   0.000  1.00  0.00           C
ATOM      8  O   GLY A   2       4.770  -1.220   0.000  1.00  0.00           O
ATOM      9  N   GLY A   3       6.050   0.620   0.000  1.00  0.00           N
ATOM     10  CA  GLY A   3       7.310  -0.080   0.000  1.00  0.00           C
ATOM     11  C   GLY A   3       8.500   0.870   0.000  1.00  0.00           C
ATOM     12  O   GLY A   3       8.380   2.090   0.000  1.00  0.00           O
END
";

impl Default for BiostructPanel {
    fn default() -> Self {
        BiostructPanel {
            tool: Tool::Analyze,
            structure_a: DEMO_PDB.to_string(),
            structure_b: DEMO_PDB.to_string(),
            clash_tolerance: 0.4,
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
            representation: Representation::default(),
            color_scheme: ColorScheme::default(),
            molviz_params: MolvizParams::default(),
        }
    }
}

/// Collect every Cα atom coordinate from the first model of a structure.
fn ca_coords(s: &Structure) -> Vec<Point3<f64>> {
    let mut out = Vec::new();
    for chain in &s.first_model().chains {
        for res in &chain.residues {
            if let Some(ca) = res.ca() {
                out.push(ca.coord);
            }
        }
    }
    out
}

/// Render the Macromolecular Structure panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.biostruct;

    common::section(ui, "Tool");
    ui.horizontal_wrapped(|ui| {
        ui.selectable_value(&mut p.tool, Tool::Analyze, "Structure analysis")
            .on_hover_text("Detect secondary structure, contacts, clashes, and chains.");
        ui.selectable_value(&mut p.tool, Tool::Ramachandran, "Ramachandran")
            .on_hover_text(
                "Compute φ/ψ backbone dihedrals and classify into Ramachandran regions.",
            );
        ui.selectable_value(&mut p.tool, Tool::Superpose, "RMSD / superpose")
            .on_hover_text("Superpose two structures via Kabsch rotation + report RMSD.");
        ui.separator();
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u {
            p.undo_edit();
        }
        if r {
            p.redo_edit();
        }
    });
    ui.separator();

    match p.tool {
        Tool::Analyze | Tool::Ramachandran => draw_single(p, ui),
        Tool::Superpose => draw_superpose(p, ui),
    }

    common::error_line(ui, &p.error);

    // --- 3-D viewport integration ---------------------------------
    // The structure-A text feeds the viewer for every tool (it is the
    // analyse / Ramachandran input and the mobile superpose
    // structure). A reactive representation picker selects how it is
    // meshed (ball-and-stick / sticks / spacefill / cartoon / surface),
    // and "Show in 3D viewport" pushes that mesh into the app's wgpu 3-D
    // viewport. The picker is a row of named `selectable_value` widgets,
    // so an agent can switch representation through the accessibility
    // tree by widget name.
    if !app.genetics.biostruct.structure_a.trim().is_empty() {
        draw_representation_picker(app, ui);
        if ui.button("Show in 3D viewport").clicked() {
            show_in_viewport(app);
        }
    }

    let p = &app.genetics.biostruct;
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "biostruct_result", &p.result, 16);
    }
}

/// Draw the reactive **representation picker** for the 3-D viewer: a row of
/// named `selectable_value` buttons (one per [`Representation`]) plus the
/// per-representation tuning controls. The named widgets are what makes the
/// representation AI-drivable — an agent flips representation by selecting the
/// button by its `label()` through the accessibility tree.
fn draw_representation_picker(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.biostruct;
    common::section(ui, "3-D representation");
    ui.horizontal_wrapped(|ui| {
        for rep in Representation::ALL {
            ui.selectable_value(&mut p.representation, rep, rep.label())
                .on_hover_text(match rep {
                    Representation::BallAndStick => "Element spheres + bond cylinders.",
                    Representation::Sticks => "Bonds only (licorice).",
                    Representation::Spacefill => "Full van-der-Waals spheres (CPK).",
                    Representation::Cartoon => {
                        "Smooth Catmull-Rom round tube through the Cα backbone \
                         (helices/strands fatten via DSSP)."
                    }
                    Representation::Ribbon => {
                        "Flat elliptical ribbon swept along the Cα spline — a wide, \
                         thin band (helices/strands widen via DSSP)."
                    }
                    Representation::Surface => "Marching-cubes molecular surface (union-of-balls).",
                    Representation::Density => {
                        "Marching-cubes isosurface of a Gaussian electron-density-like \
                         field (sum of per-atom Gaussians; a smooth-blob model, not real QM)."
                    }
                });
        }
    });

    // Colour-scheme picker: a row of named `selectable_value` buttons (one per
    // [`ColorScheme`]). A non-default scheme (chain / residue / B-factor) makes
    // "Show in 3D viewport" build a per-vertex-coloured mesh and upload the
    // colours so the viewport renders the scheme. Each button is `labelled_by`
    // the "Colour scheme" caption so it carries an unambiguous accessibility /
    // UI-Automation Name (AI-drivable: an agent flips the scheme by selecting
    // the button by its `label()`), and the caption itself names the group.
    ui.horizontal_wrapped(|ui| {
        let caption = ui.label("Colour scheme");
        for scheme in ColorScheme::ALL {
            ui.selectable_value(&mut p.color_scheme, scheme, scheme.label())
                .labelled_by(caption.id)
                .on_hover_text(match scheme {
                    ColorScheme::Element => "CPK by element (C grey, N blue, O red, S yellow, …).",
                    ColorScheme::Chain => "A distinct hue per chain.",
                    ColorScheme::Residue => "Rainbow ramp by residue index (N→C terminus).",
                    ColorScheme::BFactor => "Blue→white→red ramp by B-factor (low→high).",
                });
        }
    });
    if p.color_scheme != ColorScheme::Element {
        // Sticks / cartoon / ribbon / surface / density have no per-atom colour
        // builder, so they take a single scheme-derived tint rather than
        // per-atom colour. Tell the user which they'll get.
        let per_atom = matches!(
            p.representation,
            Representation::BallAndStick | Representation::Spacefill
        );
        if !per_atom {
            ui.label(
                egui::RichText::new(format!(
                    "note: {} renders a single {}-derived tint (per-atom colour needs \
                     ball-and-stick or spacefill)",
                    p.representation.label().to_ascii_lowercase(),
                    p.color_scheme.label().to_ascii_lowercase(),
                ))
                .weak()
                .italics(),
            );
        }
    }

    // Per-representation tuning, only shown when relevant.
    match p.representation {
        Representation::Surface => {
            ui.horizontal(|ui| {
                ui.label("Surface grid:");
                ui.add(egui::Slider::new(&mut p.molviz_params.grid_max, 16..=128).text("cells"))
                    .on_hover_text(
                        "Marching-cubes resolution along the longest axis — higher is \
                     smoother but costs O(n³).",
                    );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.probe_radius)
                        .speed(0.05)
                        .range(0.0..=3.0)
                        .prefix("probe Å "),
                )
                .on_hover_text("Probe radius inflating each atom before the isosurface.");
            });
        }
        Representation::Density => {
            ui.horizontal(|ui| {
                ui.label("Density grid:");
                ui.add(
                    egui::Slider::new(&mut p.molviz_params.density_grid_max, 16..=128)
                        .text("cells"),
                )
                .on_hover_text(
                    "Marching-cubes resolution along the longest axis — higher is \
                     smoother but costs O(n³).",
                );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.density_sigma)
                        .speed(0.05)
                        .range(0.2..=4.0)
                        .prefix("σ Å "),
                )
                .on_hover_text(
                    "Gaussian width per atom — larger σ gives fatter, more-merged blobs.",
                );
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.density_iso)
                        .speed(0.02)
                        .range(0.05..=0.95)
                        .prefix("iso "),
                )
                .on_hover_text(
                    "Iso-level as a fraction of one atom's peak — lower wraps more of \
                     each Gaussian tail. (Phenomenological blob, not real electron density.)",
                );
            });
            ui.checkbox(
                &mut p.molviz_params.density_weight_by_element,
                "Weight density by element",
            )
            .on_hover_text(
                "Scale each atom's Gaussian by a crude electron count so heavy atoms read denser.",
            );
        }
        Representation::Cartoon => {
            ui.horizontal(|ui| {
                ui.label("Tube radius (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.tube_radius)
                        .speed(0.02)
                        .range(0.1..=2.0),
                );
            });
        }
        Representation::Ribbon => {
            ui.horizontal(|ui| {
                ui.label("Ribbon width (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.ribbon_width)
                        .speed(0.05)
                        .range(0.2..=4.0),
                )
                .on_hover_text("Half-width of the flat band (the wide axis).");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.ribbon_thickness)
                        .speed(0.02)
                        .range(0.05..=1.0)
                        .prefix("thick Å "),
                )
                .on_hover_text("Half-thickness of the band (the thin axis).");
            });
        }
        Representation::BallAndStick | Representation::Sticks => {
            ui.horizontal(|ui| {
                ui.label("Stick radius (Å):");
                ui.add(
                    egui::DragValue::new(&mut p.molviz_params.stick_radius)
                        .speed(0.01)
                        .range(0.02..=0.6),
                );
            });
        }
        Representation::Spacefill => {}
    }
}

/// Build the structure-A text into the selected [`Representation`]'s mesh and
/// push it into the app's 3-D viewport. For the cartoon, the Cα backbone and a
/// per-residue DSSP secondary-structure track are extracted first; for the
/// other modes the atom/bond [`ViewMolecule`] is meshed.
fn show_in_viewport(app: &mut ValenxApp) {
    let rep = app.genetics.biostruct.representation;
    let scheme = app.genetics.biostruct.color_scheme;
    let params = app.genetics.biostruct.molviz_params;
    match read_structure(&app.genetics.biostruct.structure_a, "viewer") {
        Ok(s) => {
            let backbone = if rep.needs_backbone() {
                ca_backbone(&s)
            } else {
                Vec::new()
            };
            // A cartoon / ribbon needs a real Cα trace — fail loudly rather
            // than silently showing nothing.
            if rep.needs_backbone() && backbone.len() < 2 {
                app.genetics.biostruct.error = Some(format!(
                    "{} needs a protein backbone (≥ 2 Cα atoms) — this \
                     structure has none; try ball-and-stick or surface",
                    rep.label().to_ascii_lowercase()
                ));
                return;
            }
            let view = ViewMolecule::from_biostruct(&s);
            let label = format!("structure.{}.{}", rep.token(), scheme.token());
            // Default (Element / CPK) keeps the original monochrome-metal path
            // (`show_molecule`). A non-default scheme builds the colour-aware
            // mesh and uploads per-vertex colours so the viewport renders the
            // scheme instead of a single material.
            let result = if scheme == ColorScheme::Element {
                let mesh = molviz::build_mesh(&view, rep, &backbone, &params);
                molecule_view::show_molecule(app, mesh, &label)
            } else {
                // Per-atom annotations (chain / residue / B-factor) in lockstep
                // with `view.atoms`, read from the structure the viewer parsed.
                let attrs = structure_atom_attrs(&s);
                let (mesh, per_tri_colors) =
                    molviz::build_mesh_colored(&view, rep, &backbone, &params, scheme, &attrs);
                molecule_view::show_molecule_colored(app, mesh, &per_tri_colors, &label)
            };
            match result {
                Ok(_) => app.genetics.biostruct.error = None,
                Err(e) => app.genetics.biostruct.error = Some(e),
            }
        }
        Err(e) => app.genetics.biostruct.error = Some(e.to_string()),
    }
}

/// Build the per-atom [`AtomAttr`] (chain id, residue index, B-factor) for a
/// structure's first model, **in the exact atom order**
/// [`ViewMolecule::from_biostruct`] walks (chain-major, then residue, then atom)
/// so the returned vec is in lockstep with `ViewMolecule::atoms` and the colour
/// schemes line up atom-for-atom.
///
/// All three fields come straight from the `valenx-biostruct` model, which
/// carries them on every atom (`Atom::b_factor`, `Chain::id`, `Residue::seq_num`):
/// no field needs a fallback. The residue index is a **monotone counter across
/// the whole model** (incremented per residue in iteration order) rather than the
/// raw `seq_num`, so the residue rainbow runs cleanly N→C even when `seq_num`
/// has gaps, insertion codes, or resets between chains.
fn structure_atom_attrs(s: &Structure) -> Vec<AtomAttr> {
    let mut out = Vec::new();
    let mut residue_index: i32 = 0;
    for chain in &s.first_model().chains {
        for res in &chain.residues {
            for atom in &res.atoms {
                out.push(AtomAttr::new(
                    chain.id.clone(),
                    residue_index,
                    atom.b_factor as f32,
                ));
            }
            residue_index += 1;
        }
    }
    out
}

/// Extract the per-residue Cα backbone control points of a structure's first
/// model, tagged with their DSSP secondary-structure code, for the cartoon
/// representation. Residues without a Cα (ligands, waters, missing atoms) are
/// skipped; the DSSP track is computed per chain via [`valenx_biostruct::dssp`]
/// and indexed by residue position so each kept Cα carries its own state.
fn ca_backbone(s: &valenx_biostruct::structure::Structure) -> Vec<BackbonePoint> {
    use valenx_biostruct::dssp;
    let mut out = Vec::new();
    for chain in &s.first_model().chains {
        // Per-chain DSSP state, one entry per residue in chain order.
        let states = dssp::assign_chain(chain).states;
        for (i, res) in chain.residues.iter().enumerate() {
            if let Some(ca) = res.ca() {
                let ss = states.get(i).map(|st| st.code());
                out.push(BackbonePoint::new(
                    [ca.coord.x as f32, ca.coord.y as f32, ca.coord.z as f32],
                    ss,
                ));
            }
        }
    }
    out
}

fn structure_text_input(
    ui: &mut egui::Ui,
    id: &str,
    label: &str,
    buf: &mut String,
) -> Option<String> {
    let mut err = None;
    common::section(ui, label);
    ui.horizontal(|ui| {
        if ui.small_button("Load PDB / mmCIF…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Structure", &["pdb", "ent", "cif", "mmcif"])
                .pick_file()
            {
                // Round-21 H1: file-dialog paths flow straight to a
                // bare `fs::read_to_string` pre-fix. A user (or a
                // stale dialog state) pointing at a multi-GB file
                // would OOM the renderer before the parser saw a
                // single byte. `read_capped_to_string` rejects
                // anything past `MAX_GENETICS_FILE_BYTES` (64 MiB).
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => *buf = t,
                    Err(e) => err = Some(format!("read: {e}")),
                }
            }
        }
        ui.label(format!("{} lines", buf.lines().count()));
    });
    ui.add(
        egui::TextEdit::multiline(buf)
            .id_source(id)
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(6),
    );
    err
}

fn draw_single(p: &mut BiostructPanel, ui: &mut egui::Ui) {
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_a",
        "Structure (PDB / mmCIF)",
        &mut p.structure_a,
    ) {
        p.error = Some(e);
    }
    if p.tool == Tool::Analyze {
        ui.horizontal(|ui| {
            ui.label("Clash tolerance (Å):");
            ui.add(
                egui::DragValue::new(&mut p.clash_tolerance)
                    .speed(0.05)
                    .range(0.0..=2.0),
            );
        });
    }
    let label = if p.tool == Tool::Analyze {
        "Analyze structure"
    } else {
        "Classify Ramachandran"
    };
    if common::run_button(ui, label) {
        let snap = p.snapshot();
        p.history.record(snap);
        run_single(p);
    }
}

/// Run the structure-analysis / Ramachandran tool — extracted from the
/// button closure so it is callable from the headless UI tests.
fn run_single(p: &mut BiostructPanel) {
    p.error = None;
    match read_structure(&p.structure_a, "input") {
        Ok(s) => match p.tool {
            Tool::Analyze => match StructureReport::analyze(&s, p.clash_tolerance) {
                Ok(r) => {
                    let mut out = format!(
                        "title          : {}\nmodels         : {}\n\
                             chains         : {}  ({} protein, {} nucleic)\n\
                             residues       : {}\natoms          : {}\n\
                             water / ligand : {} / {}\nradius of gyr. : {:.2} Å\n\
                             mean helix     : {:.1} %\nmean sheet     : {:.1} %\n\n\
                             -- per chain --\n",
                        r.title,
                        r.model_count,
                        r.chains.len(),
                        r.protein_chain_count(),
                        r.nucleic_chain_count(),
                        r.residue_count,
                        r.atom_count,
                        r.water_count,
                        r.ligand_count,
                        r.radius_of_gyration,
                        r.mean_helix_fraction() * 100.0,
                        r.mean_sheet_fraction() * 100.0,
                    );
                    for ch in &r.chains {
                        out.push_str(&format!(
                            "  {} {:<10} {:>4} res  H {:>4.0}% E {:>4.0}% C {:>4.0}%\n",
                            ch.id,
                            format!("{:?}", ch.kind),
                            ch.residue_count,
                            ch.secondary.helix * 100.0,
                            ch.secondary.sheet * 100.0,
                            ch.secondary.coil * 100.0,
                        ));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            },
            Tool::Ramachandran => {
                let mut out = String::new();
                for chain in &s.first_model().chains {
                    let summary = rama_summarize(chain);
                    if summary.total == 0 {
                        continue;
                    }
                    out.push_str(&format!(
                        "chain {} — {} phi/psi points\n  alpha-helix : {}\n  \
                             beta-sheet  : {}\n  left-alpha  : {}\n  bridge      : {}\n  \
                             outliers    : {}\n  allowed     : {:.1} %\n\n",
                        chain.id,
                        summary.total,
                        summary.alpha,
                        summary.beta,
                        summary.left_alpha,
                        summary.bridge,
                        summary.outliers,
                        summary.allowed_fraction() * 100.0,
                    ));
                }
                if out.is_empty() {
                    out = "no residues with a defined φ/ψ pair (need ≥ 3 \
                               consecutive amino acids)"
                        .to_string();
                }
                p.result = out;
            }
            Tool::Superpose => unreachable!(),
        },
        Err(e) => p.error = Some(e.to_string()),
    }
}

fn draw_superpose(p: &mut BiostructPanel, ui: &mut egui::Ui) {
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_mob",
        "Mobile structure",
        &mut p.structure_a,
    ) {
        p.error = Some(e);
    }
    if let Some(e) = structure_text_input(
        ui,
        "biostruct_input_ref",
        "Reference structure",
        &mut p.structure_b,
    ) {
        p.error = Some(e);
    }
    if common::run_button(ui, "Kabsch superpose (Cα)") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_superpose(p);
    }
}

/// Run the Kabsch Cα superposition — extracted for the headless UI
/// tests.
fn run_superpose(p: &mut BiostructPanel) {
    p.error = None;
    match (
        read_structure(&p.structure_a, "mobile"),
        read_structure(&p.structure_b, "reference"),
    ) {
        (Ok(mob), Ok(reference)) => {
            let ca_m = ca_coords(&mob);
            let ca_r = ca_coords(&reference);
            let n = ca_m.len().min(ca_r.len());
            if n < 3 {
                p.error = Some(format!(
                    "need ≥ 3 paired Cα atoms (mobile {}, reference {})",
                    ca_m.len(),
                    ca_r.len(),
                ));
                return;
            }
            let m = &ca_m[..n];
            let r = &ca_r[..n];
            let pre = rmsd(m, r);
            match kabsch(m, r) {
                Ok(sup) => {
                    p.result = format!(
                        "paired Cα atoms : {}\nRMSD before     : {}\n\
                         RMSD after fit  : {:.4} Å\n\n\
                         optimal rotation + translation found by the \
                         Kabsch algorithm.",
                        n,
                        pre.map(|v| format!("{v:.4} Å"))
                            .unwrap_or_else(|_| "(n/a)".into()),
                        sup.rmsd,
                    );
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        (Err(e), _) => p.error = Some(format!("mobile: {e}")),
        (_, Err(e)) => p.error = Some(format!("reference: {e}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demo_pdb_parses() {
        let s = read_structure(DEMO_PDB, "demo").expect("demo PDB must parse");
        assert_eq!(ca_coords(&s).len(), 3);
    }

    #[test]
    fn demo_superposes_to_zero_rmsd() {
        // The same structure superposed on itself has ~zero RMSD.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let ca = ca_coords(&s);
        let sup = kabsch(&ca, &ca).unwrap();
        assert!(sup.rmsd < 1.0e-6);
    }

    /// Round-21 H1 RED→GREEN: the genetics-workbench file loaders
    /// route through [`valenx_core::io_caps::read_capped_to_string`]
    /// with the [`valenx_core::io_caps::MAX_GENETICS_FILE_BYTES`]
    /// cap. Pre-fix the loader did a bare `fs::read_to_string`, so a
    /// user-picked multi-GB file would OOM. We exercise the helper
    /// here against an oversized scratch file (allocating 100 MiB on
    /// disk would slow CI, so the test uses a small cap so a
    /// modestly-sized file trips it the same way 64 MiB+ would trip
    /// the production cap).
    #[test]
    fn oversize_genetics_load_returns_invalid_data() {
        use std::io::Write;
        let tmp = std::env::temp_dir().join("valenx_r21_biostruct_oversize.pdb");
        let mut f = std::fs::File::create(&tmp).unwrap();
        f.write_all(&vec![b'A'; 4096]).unwrap();
        drop(f);
        // Cap of 1 KiB simulates the 64 MiB production cap shape —
        // a file larger than the cap is rejected with InvalidData.
        let err = valenx_core::io_caps::read_capped_to_string(&tmp, 1024).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
        let _ = std::fs::remove_file(&tmp);
        // Sanity: the constant the loader actually uses is the
        // production cap (proves we didn't downgrade by accident).
        assert_eq!(
            valenx_core::io_caps::MAX_GENETICS_FILE_BYTES,
            64u64 * 1024 * 1024
        );
    }
}

/// Headless egui UI-logic tests for the Macromolecular Structure panel.
#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use crate::genetics_workbench::GeneticsPanel;
    use crate::ValenxApp;

    fn draw_headless(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
    }

    /// Draw the panel once with the accesskit tree enabled and return its
    /// nodes — the harness the AI-drivability (`labelled_by`) assertions use.
    fn draw_and_collect_nodes(
        app: &mut ValenxApp,
    ) -> Vec<(egui::accesskit::NodeId, egui::accesskit::Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(app, ui);
            });
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::MacromolecularStructure;
        app
    }

    #[test]
    fn draws_every_tool_without_panic() {
        for tool in [Tool::Analyze, Tool::Ramachandran, Tool::Superpose] {
            let mut app = app_with_panel();
            app.genetics.biostruct.tool = tool;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.biostruct.result = "chains : 1\nresidues : 3\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.biostruct.error = Some("could not parse PDB".to_string());
        draw_headless(&mut app);
        // Empty structure text — the "Show in 3D" affordance is gated.
        let mut app = app_with_panel();
        app.genetics.biostruct.structure_a.clear();
        draw_headless(&mut app);
    }

    #[test]
    fn run_single_analyzes_the_demo_structure() {
        // The demo glycine peptide → the real valenx-biostruct
        // analyzer produces a correctly-formatted report.
        let mut p = BiostructPanel {
            tool: Tool::Analyze,
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_none(), "analyze errored: {:?}", p.error);
        assert!(p.result.contains("chains"));
        assert!(p.result.contains("residues"));
    }

    #[test]
    fn run_single_classifies_ramachandran() {
        let mut p = BiostructPanel {
            tool: Tool::Ramachandran,
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_none(), "Ramachandran errored: {:?}", p.error);
        assert!(!p.result.is_empty());
    }

    #[test]
    fn run_superpose_aligns_identical_structures() {
        // The demo structure superposed on itself → near-zero RMSD.
        let mut p = BiostructPanel {
            tool: Tool::Superpose,
            ..BiostructPanel::default()
        };
        run_superpose(&mut p);
        assert!(p.error.is_none(), "superpose errored: {:?}", p.error);
        assert!(p.result.contains("paired Cα atoms"));
        assert!(p.result.contains("RMSD after fit"));
    }

    #[test]
    fn run_actions_surface_errors_on_bad_input() {
        // A non-PDB string is malformed structure input.
        let mut p = BiostructPanel {
            tool: Tool::Analyze,
            structure_a: "not a structure".to_string(),
            ..BiostructPanel::default()
        };
        run_single(&mut p);
        assert!(p.error.is_some(), "analyze should error on malformed input");
        // Superpose with too few Cα atoms (an empty mobile structure).
        let mut p = BiostructPanel {
            tool: Tool::Superpose,
            structure_a: "END\n".to_string(),
            ..BiostructPanel::default()
        };
        run_superpose(&mut p);
        assert!(p.error.is_some(), "superpose should error with no Cα atoms");
    }

    #[test]
    fn draws_with_every_representation_selected_without_panic() {
        // The representation picker + its per-mode tuning row must draw for
        // every representation (the demo PDB is loaded by default, so the
        // picker is shown).
        for rep in Representation::ALL {
            let mut app = app_with_panel();
            app.genetics.biostruct.representation = rep;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn show_in_viewport_meshes_the_demo_for_each_representation() {
        // The demo glycine peptide is a real protein backbone, so every
        // representation — including cartoon (needs Cα) and surface (marching
        // cubes) — produces a non-empty mesh in the viewport without error.
        for rep in Representation::ALL {
            let mut app = app_with_panel();
            app.genetics.biostruct.representation = rep;
            // Keep the surface grid small so the test stays fast.
            app.genetics.biostruct.molviz_params.grid_max = 24;
            super::show_in_viewport(&mut app);
            assert!(
                app.genetics.biostruct.error.is_none(),
                "{rep:?} errored: {:?}",
                app.genetics.biostruct.error
            );
            assert!(
                app.stl.is_some(),
                "{rep:?} should have pushed a mesh to the viewport"
            );
        }
    }

    #[test]
    fn cartoon_errors_on_a_structure_with_no_backbone() {
        // A hetero-only structure (a lone zinc ion) has no Cα trace, so the
        // cartoon representation must fail loudly rather than show nothing.
        let mut app = app_with_panel();
        app.genetics.biostruct.representation = Representation::Cartoon;
        app.genetics.biostruct.structure_a = "\
HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00  0.00          ZN
END
"
        .to_string();
        super::show_in_viewport(&mut app);
        assert!(
            app.genetics.biostruct.error.is_some(),
            "cartoon with no backbone should surface an error"
        );
    }

    #[test]
    fn ca_backbone_extracts_dssp_tagged_trace() {
        // The demo 3-glycine peptide yields 3 Cα control points, each carrying
        // a DSSP secondary-structure code (so the cartoon tube can modulate).
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let bb = super::ca_backbone(&s);
        assert_eq!(bb.len(), 3, "three Cα control points");
        assert!(
            bb.iter().all(|p| p.ss.is_some()),
            "every backbone point carries a DSSP code"
        );
    }

    #[test]
    fn draws_with_every_color_scheme_selected_without_panic() {
        // The colour-scheme picker (+ its per-scheme note row) must draw for
        // every scheme, in combination with each representation (the demo PDB
        // is loaded by default, so the picker is shown).
        for scheme in ColorScheme::ALL {
            for rep in Representation::ALL {
                let mut app = app_with_panel();
                app.genetics.biostruct.color_scheme = scheme;
                app.genetics.biostruct.representation = rep;
                draw_headless(&mut app);
            }
        }
    }

    #[test]
    fn structure_atom_attrs_lockstep_with_view_atoms() {
        // The per-atom attribute slice must be in lockstep with the
        // `ViewMolecule::from_biostruct` atom order so the colour schemes line
        // up atom-for-atom. The demo 3-glycine peptide has 12 atoms across one
        // chain "A" and three residues.
        let s = read_structure(DEMO_PDB, "demo").unwrap();
        let view = ViewMolecule::from_biostruct(&s);
        let attrs = super::structure_atom_attrs(&s);
        assert_eq!(
            attrs.len(),
            view.atoms.len(),
            "one AtomAttr per ViewMolecule atom, in order"
        );
        assert!(attrs.iter().all(|a| a.chain == "A"), "demo is all chain A");
        // Residue index is a monotone 0-based counter across the model: three
        // residues → indices 0, 1, 2 present.
        let max_res = attrs.iter().map(|a| a.residue_index).max().unwrap();
        assert_eq!(max_res, 2, "three residues → max residue index 2");
    }

    #[test]
    fn show_in_viewport_uploads_per_vertex_colors_for_each_non_default_scheme() {
        // A non-default colour scheme builds a colour-aware mesh and uploads
        // per-vertex colours: the viewport STL carries a `colors` buffer whose
        // length is exactly 3 × the triangle count (one colour per surface
        // vertex), for every representation. The default Element scheme keeps
        // the monochrome path (no colour buffer).
        for scheme in ColorScheme::ALL {
            for rep in Representation::ALL {
                let mut app = app_with_panel();
                app.genetics.biostruct.color_scheme = scheme;
                app.genetics.biostruct.representation = rep;
                // Keep the surface/density grids small so the test stays fast.
                app.genetics.biostruct.molviz_params.grid_max = 24;
                app.genetics.biostruct.molviz_params.density_grid_max = 24;
                super::show_in_viewport(&mut app);
                assert!(
                    app.genetics.biostruct.error.is_none(),
                    "{scheme:?}/{rep:?} errored: {:?}",
                    app.genetics.biostruct.error
                );
                let stl = app
                    .stl
                    .as_ref()
                    .unwrap_or_else(|| panic!("{scheme:?}/{rep:?} pushed no mesh"));
                let tri_count = stl.mesh.triangle_count();
                match scheme {
                    ColorScheme::Element => {
                        // Default scheme → no per-vertex colour override.
                        assert!(
                            stl.colors.is_none(),
                            "{rep:?}: Element scheme must keep the monochrome path"
                        );
                    }
                    _ => {
                        let colors = stl.colors.as_ref().unwrap_or_else(|| {
                            panic!("{scheme:?}/{rep:?} must carry per-vertex colours")
                        });
                        assert_eq!(
                            colors.len(),
                            tri_count * 3,
                            "{scheme:?}/{rep:?}: colours must equal 3 × triangle count"
                        );
                        // Every colour component is a finite 0..=1 value.
                        assert!(
                            colors.iter().all(|c| c
                                .iter()
                                .all(|&x| x.is_finite() && (0.0..=1.0).contains(&x))),
                            "{scheme:?}/{rep:?}: colour components must be finite in [0,1]"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn color_scheme_picker_is_labelled_for_ai_driving() {
        // The colour-scheme picker must be AI-drivable: each scheme button
        // carries the scheme label as its accessibility Name AND is
        // `labelled_by` the "Colour scheme" caption (so an agent / screen
        // reader can find and select it by name). Assert the four scheme labels
        // appear as named, `labelled_by`-associated nodes in the accesskit tree.
        use egui::accesskit::Node;
        let mut app = app_with_panel();
        let nodes = draw_and_collect_nodes(&mut app);

        for scheme in ColorScheme::ALL {
            let label = scheme.label();
            let matching: Vec<&Node> = nodes
                .iter()
                .map(|(_, n)| n)
                .filter(|n| n.name() == Some(label))
                .collect();
            assert!(
                !matching.is_empty(),
                "colour scheme '{label}' must appear as a named node in the accesskit tree"
            );
            assert!(
                matching.iter().any(|n| !n.labelled_by().is_empty()),
                "colour scheme '{label}' button must be labelled_by its caption (AI-drivable)"
            );
        }
    }
}
