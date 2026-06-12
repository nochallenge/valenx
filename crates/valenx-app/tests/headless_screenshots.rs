//! Headless screenshot harness — render every workbench panel +
//! visualization mode to PNG files so a developer / reviewer has
//! concrete artifacts to inspect for the "live visual aesthetic"
//! check without running the app on a machine.
//!
//! # What this proves
//!
//! Valenx already has 151+ headless UI-logic tests (the
//! `headless_ui_tests` name-filtered set) plus a verified headless GPU
//! render path (the PBR shader render-with-readback test in
//! `wgpu_renderer.rs`). What was missing: nobody could *see* what the
//! workbench panels actually look like without running the binary.
//!
//! This integration test closes that gap. For each of ~33 workbench
//! panels (Mesh/CAD toolbox sub-panels, Genetics workbench panels, Aero
//! / Wind Tunnel workflow sections) it:
//!
//! 1. Builds a `ValenxApp` with the panel toggled on and any state the
//!    populated-state populators expect.
//! 2. Runs one egui frame with the panel's draw function inside a
//!    `CentralPanel` (or, for the host SidePanels, the whole workbench).
//! 3. Tessellates the shapes, hands the paint jobs to a real
//!    `egui_wgpu::Renderer` pointed at an off-screen colour texture.
//! 4. Encodes the render pass, submits, copies the texture to a 256-
//!    byte-row-aligned buffer, maps it, and writes the pixels out as a
//!    PNG to `screenshots/<workbench>/<panel>.png`.
//!
//! # Honest environment handling
//!
//! `wgpu::Instance::request_adapter` returns `None` in a sandbox with
//! no GPU and no software fallback. We log "no adapter, skipping" and
//! return cleanly — the test passes (no PNGs written, no crash).
//!
//! # What this does NOT do
//!
//! - It does not perform an aesthetic judgement. That's the human's job
//!   reviewing the PNGs.
//! - It does not exercise every interaction state, every theme, every
//!   locale. One representative frame per panel.
//! - It does not screenshot the central 3-D viewport or the wgpu PBR
//!   pass — those are validated by the existing headless render tests
//!   in `wgpu_renderer.rs` (`headless_pbr_render_shades_a_lit_quad`).
//!
//! # Running
//!
//! ```text
//! cargo test -p valenx-app --test headless_screenshots
//! ```
//!
//! Outputs land in `screenshots/` at the workspace root. The directory
//! is `.gitignore`-excluded — the harness regenerates them on demand.

// `push_cad_sub_panel` / `push_aero_section` take a `sub_draw` fn
// pointer as load-bearing in-line documentation right next to the
// catalogue's slug — the actual dispatch goes through `select_*_draw`
// by slug, so clippy flags the duplicate `|a, u| f(a, u)` adapters as
// "redundant closure". We intentionally keep them for the visual
// catalogue + the type-check.
#![allow(clippy::redundant_closure)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use eframe::egui;
use eframe::egui_wgpu;
use eframe::wgpu;

use valenx_app::aero_workbench::draw_aero_workbench;
use valenx_app::genetics_workbench::{draw_genetics_workbench, GeneticsPanel};
use valenx_app::mesh_toolbox::draw_mesh_toolbox;
use valenx_app::ValenxApp;

// ---------------------------------------------------------------------
// Test entry point
// ---------------------------------------------------------------------

/// Drive the full screenshot harness. Builds an off-screen wgpu
/// device, walks the panel catalogue, and writes one PNG per panel.
///
/// Skips cleanly (no failure, no PNGs) when no wgpu adapter exists.
#[test]
fn render_every_workbench_panel_to_png() {
    let Some(mut harness) = ScreenshotHarness::new() else {
        eprintln!(
            "render_every_workbench_panel_to_png: no wgpu adapter — \
             skipping (CI without GPU / no software fallback). The \
             panel headless_ui_tests still ran and proved every panel \
             draws without panicking; this harness only adds the PNG \
             artifact."
        );
        return;
    };

    let out_root = workspace_root().join("screenshots");
    let _ = fs::remove_dir_all(&out_root);
    fs::create_dir_all(&out_root).expect("create screenshots/");

    let mut all_renders: Vec<RenderedPanel> = Vec::new();
    for panel in panel_catalogue() {
        let dir = out_root.join(panel.workbench_slug());
        fs::create_dir_all(&dir).expect("create workbench dir");
        let out = dir.join(format!("{}.png", panel.slug));

        let dims = harness.render_panel_to_png(&panel, &out);

        all_renders.push(RenderedPanel {
            spec: panel,
            out_path: out,
            size: dims,
        });
    }

    // Per-file assertions: every PNG should exist, be > 1 KB, have the
    // right dimensions, and not be uniformly the clear colour.
    for rendered in &all_renders {
        let meta = fs::metadata(&rendered.out_path).unwrap_or_else(|e| {
            panic!(
                "screenshot for `{}` was never written ({}): {e}",
                rendered.spec.slug,
                rendered.out_path.display(),
            )
        });
        assert!(
            meta.len() > 1024,
            "screenshot for `{}` is suspiciously small ({} bytes) — \
             rendering likely produced an empty / uniform frame: {}",
            rendered.spec.slug,
            meta.len(),
            rendered.out_path.display(),
        );
        assert_eq!(
            rendered.size,
            [PANEL_WIDTH, PANEL_HEIGHT],
            "screenshot for `{}` rendered at unexpected dimensions",
            rendered.spec.slug,
        );
    }

    write_index(&out_root, &all_renders).expect("write INDEX.md");

    eprintln!(
        "render_every_workbench_panel_to_png: wrote {} PNGs to {}",
        all_renders.len(),
        out_root.display()
    );
}

// ---------------------------------------------------------------------
// Texture dimensions — kept reasonable so the harness finishes quickly
// (~1 MB per PNG raw, ~100 KB encoded) and the dimensions assertion
// holds.
// ---------------------------------------------------------------------

const PANEL_WIDTH: u32 = 1280;
const PANEL_HEIGHT: u32 = 800;

// ---------------------------------------------------------------------
// Panel catalogue — every workbench panel we render, with a slug, a
// human-readable label, and a one-line description for the INDEX.
// ---------------------------------------------------------------------

/// What kind of root container the panel needs.
#[derive(Copy, Clone, Debug)]
enum Mount {
    /// Run the whole host SidePanel (e.g. the entire Mesh Toolbox right
    /// panel). The egui frame mounts via the host's `draw_*` fn against
    /// `egui::Context` directly.
    HostSidePanel,
    /// Run an individual panel inside a `CentralPanel` (the same pattern
    /// every `headless_ui_tests` `draw_headless` uses).
    CentralPanel,
}

/// What workbench the panel belongs to. Drives the output sub-directory.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum Workbench {
    MeshCad,
    Genetics,
    Aero,
}

impl Workbench {
    fn slug(self) -> &'static str {
        match self {
            Workbench::MeshCad => "mesh-cad",
            Workbench::Genetics => "genetics",
            Workbench::Aero => "aero",
        }
    }
    fn label(self) -> &'static str {
        match self {
            Workbench::MeshCad => "Mesh / CAD Toolbox",
            Workbench::Genetics => "Genetics Workbench",
            Workbench::Aero => "Aerodynamics / Wind Tunnel",
        }
    }
}

/// One panel-render request — what to draw, where to put the PNG.
struct PanelSpec {
    workbench: Workbench,
    slug: &'static str,
    label: &'static str,
    description: &'static str,
    mount: Mount,
    /// Closure that populates a fresh `ValenxApp` to the panel's
    /// representative state and draws it into the supplied egui context.
    draw: PanelDrawFn,
}

impl PanelSpec {
    fn workbench_slug(&self) -> &'static str {
        self.workbench.slug()
    }
}

type PanelDrawFn = fn(&mut ValenxApp, &egui::Context);

struct RenderedPanel {
    spec: PanelSpec,
    out_path: PathBuf,
    size: [u32; 2],
}

/// Every panel we render — the human-curated catalogue.
fn panel_catalogue() -> Vec<PanelSpec> {
    let mut v: Vec<PanelSpec> = Vec::new();

    // ----- Mesh / CAD Toolbox ----------------------------------------
    // The whole host SidePanel with every collapsing section + every
    // sub-panel stacked inside it. The reviewer can see how the toolbox
    // looks when fully expanded.
    v.push(PanelSpec {
        workbench: Workbench::MeshCad,
        slug: "00-toolbox-host",
        label: "Mesh Toolbox (host panel)",
        description: "The whole right-side Mesh Toolbox SidePanel — inspector, \
             part workbench, transformations, cut plane, repair, export, \
             plus every CAD sub-panel as a collapsing header.",
        mount: Mount::HostSidePanel,
        draw: |app, ctx| {
            app.enable_mesh_toolbox();
            draw_mesh_toolbox(app, ctx);
        },
    });
    // Individual sub-panels inside a CentralPanel. Each driven by its
    // pub(crate) draw fn — exposed through the public module path via
    // the `headless_ui_tests` blocks under valenx-app.
    push_cad_sub_panel(
        &mut v,
        "01-part",
        "Part",
        "Part workbench — primitives (box / cylinder / sphere / cone / \
         torus) + boolean ops (union / cut / intersect) over the truck \
         BRep kernel.",
        |app, ui| valenx_app::mesh_toolbox::draw_part_design_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "02-draft",
        "Draft",
        "2-D draft workbench — line, polyline, arc, circle, rectangle, \
         polygon, dimension, text via the valenx-draft kernel.",
        |app, ui| valenx_app::mesh_toolbox::draw_draft_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "03-techdraw",
        "TechDraw",
        "Technical-drawing workbench — projection views, dimensioning, \
         leader lines, title block, hatch fills.",
        |app, ui| valenx_app::mesh_toolbox::draw_techdraw_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "04-assembly",
        "Assembly",
        "Assembly workbench — parts list, mating constraints, exploded \
         views, the valenx-assembly graph view.",
        |app, ui| valenx_app::mesh_toolbox::draw_assembly_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "05-surface",
        "Surface",
        "Free-form surface workbench — NURBS curves + surfaces, Coons \
         patch, sew, trim, knot edits, surface-surface intersection.",
        |app, ui| valenx_app::mesh_toolbox::draw_surface_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "06-cam",
        "CAM",
        "CAM workbench — stock, tool table, operations (face / contour \
         / pocket / drill / adaptive), G-code generation + post.",
        |app, ui| valenx_app::mesh_toolbox::draw_cam_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "07-arch",
        "Arch / BIM",
        "Architecture / BIM workbench — walls, floors, windows, doors, \
         IFC export.",
        |app, ui| valenx_app::mesh_toolbox::draw_arch_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "08-spreadsheet",
        "Spreadsheet",
        "Spreadsheet workbench — named cells, formulas, drive CAD \
         parameters from a spreadsheet.",
        |app, ui| valenx_app::mesh_toolbox::draw_spreadsheet_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "09-dock",
        "Dock",
        "Dock workbench — layout the active panels into split / tab \
         containers.",
        |app, ui| valenx_app::mesh_toolbox::draw_dock_panel(app, ui),
    );
    push_cad_sub_panel(
        &mut v,
        "10-sketcher",
        "Sketcher",
        "2-D sketcher workbench — constrained 2-D entities feeding the \
         Part Design workbench's profiles.",
        |app, ui| valenx_app::mesh_toolbox::draw_sketcher_panel(app, ui),
    );

    // ----- Genetics Workbench ----------------------------------------
    // The host SidePanel (defaults to the Sequence panel).
    v.push(PanelSpec {
        workbench: Workbench::Genetics,
        slug: "00-workbench-host",
        label: "Genetics Workbench (host panel)",
        description: "The whole right-side Genetics Workbench SidePanel — \
             tab selector for 14 panels, default Sequence panel shown.",
        mount: Mount::HostSidePanel,
        draw: |app, ctx| {
            app.enable_genetics_workbench(GeneticsPanel::Sequence);
            draw_genetics_workbench(app, ctx);
        },
    });
    // One PNG per genetics panel — switch `genetics.active`, draw the
    // whole workbench (the active panel renders inside).
    for (i, panel) in GeneticsPanel::ALL.iter().enumerate() {
        v.push(PanelSpec {
            workbench: Workbench::Genetics,
            slug: genetics_slug(*panel, i),
            label: genetics_label(*panel),
            description: genetics_description(*panel),
            mount: Mount::HostSidePanel,
            // We can't capture the panel by value in a fn pointer, but
            // we can dispatch through a static jump table.
            draw: genetics_draw_fn(*panel),
        });
    }

    // ----- Aero / Wind Tunnel Workbench -------------------------------
    // The host SidePanel (the eight workflow sections are all visible
    // inside one scrollable column).
    v.push(PanelSpec {
        workbench: Workbench::Aero,
        slug: "00-workbench-host",
        label: "Wind Tunnel Workbench (host panel)",
        description: "The whole right-side Wind Tunnel SidePanel — Body, Wind, \
             Ground, Tunnel, Solver, Run, Results, Visualization \
             sections all visible.",
        mount: Mount::HostSidePanel,
        draw: |app, ctx| {
            app.enable_aero_workbench();
            draw_aero_workbench(app, ctx);
        },
    });
    // Individual aero sections inside a CentralPanel.
    push_aero_section(
        &mut v,
        "01-body",
        "Body",
        "Body section — geometry source picker (built-in NACA / sphere \
         / loaded mesh), reference area + chord inputs.",
        |app, ui| valenx_app::aero::panels::draw_body_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "02-wind",
        "Wind",
        "Wind conditions section — free-stream velocity, density, \
         viscosity, angle of attack, side-slip angle.",
        |app, ui| valenx_app::aero::panels::draw_wind_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "03-ground",
        "Ground",
        "Ground & wheels section — ground plane, moving belt velocity, \
         rotating-wheel velocities.",
        |app, ui| valenx_app::aero::panels::draw_ground_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "04-tunnel",
        "Tunnel & Mesh",
        "Tunnel & mesh section — virtual tunnel dimensions, \
         immersed-boundary grid resolution, wall-spacing y+.",
        |app, ui| valenx_app::aero::panels::draw_tunnel_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "05-solver",
        "Solver",
        "Solver section — SIMPLE / PISO scheme, turbulence model \
         (k-ε / k-ω SST / laminar), iteration cap, residual tolerance.",
        |app, ui| valenx_app::aero::panels::draw_solver_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "06-run",
        "Run",
        "Run section — Start / Cancel buttons, live residual chart, \
         iteration counter, status line.",
        |app, ui| valenx_app::aero::panels::draw_run_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "07-results",
        "Results",
        "Results section — drag / lift / moment coefficients, polar \
         curve, converged residual summary.",
        |app, ui| valenx_app::aero::panels::draw_results_section(app, ui),
    );
    push_aero_section(
        &mut v,
        "08-visualization",
        "Visualization",
        "Flow visualization section — pressure / velocity / vorticity \
         field overlays, push to 3-D viewport.",
        |app, ui| valenx_app::aero::panels::draw_visualization_section(app, ui),
    );

    v
}

fn push_cad_sub_panel(
    v: &mut Vec<PanelSpec>,
    slug: &'static str,
    label: &'static str,
    description: &'static str,
    _sub_draw: fn(&mut ValenxApp, &mut egui::Ui),
) {
    // Slug → fn-pointer dispatch via `select_cad_draw`. The `_sub_draw`
    // parameter is retained at the call site so the catalogue keeps the
    // CAD sub-panel draw fns visible right next to their slugs (the
    // duplication is intentional documentation).
    v.push(PanelSpec {
        workbench: Workbench::MeshCad,
        slug,
        label,
        description,
        mount: Mount::CentralPanel,
        draw: select_cad_draw(slug),
    });
}

fn push_aero_section(
    v: &mut Vec<PanelSpec>,
    slug: &'static str,
    label: &'static str,
    description: &'static str,
    _sub_draw: fn(&mut ValenxApp, &mut egui::Ui),
) {
    v.push(PanelSpec {
        workbench: Workbench::Aero,
        slug,
        label,
        description,
        mount: Mount::CentralPanel,
        draw: select_aero_draw(slug),
    });
}

/// Look up the CAD sub-panel draw fn by slug. The catalogue above only
/// stores the slug; the actual draw fns live in this jump table.
fn select_cad_draw(slug: &'static str) -> PanelDrawFn {
    match slug {
        "01-part" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_part_design_panel(app, ui);
            })
        },
        "02-draft" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_draft_panel(app, ui);
            })
        },
        "03-techdraw" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_techdraw_panel(app, ui);
            })
        },
        "04-assembly" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_assembly_panel(app, ui);
            })
        },
        "05-surface" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_surface_panel(app, ui);
            })
        },
        "06-cam" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_cam_panel(app, ui);
            })
        },
        "07-arch" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_arch_panel(app, ui);
            })
        },
        "08-spreadsheet" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_spreadsheet_panel(app, ui);
            })
        },
        "09-dock" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_dock_panel(app, ui);
            })
        },
        "10-sketcher" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::mesh_toolbox::draw_sketcher_panel(app, ui);
            })
        },
        other => panic!("unknown CAD sub-panel slug: {other}"),
    }
}

fn select_aero_draw(slug: &'static str) -> PanelDrawFn {
    match slug {
        "01-body" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_body_section(app, ui);
            })
        },
        "02-wind" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_wind_section(app, ui);
            })
        },
        "03-ground" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_ground_section(app, ui);
            })
        },
        "04-tunnel" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_tunnel_section(app, ui);
            })
        },
        "05-solver" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_solver_section(app, ui);
            })
        },
        "06-run" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_run_section(app, ui);
            })
        },
        "07-results" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_results_section(app, ui);
            })
        },
        "08-visualization" => |app, ctx| {
            draw_in_central(ctx, |ui| {
                valenx_app::aero::panels::draw_visualization_section(app, ui);
            })
        },
        other => panic!("unknown aero section slug: {other}"),
    }
}

fn draw_in_central(ctx: &egui::Context, ui_fn: impl FnOnce(&mut egui::Ui)) {
    egui::CentralPanel::default().show(ctx, |ui| {
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, ui_fn);
    });
}

/// Dispatch the genetics panel by enum to a static fn pointer that
/// activates the right panel and draws the whole workbench.
fn genetics_draw_fn(panel: GeneticsPanel) -> PanelDrawFn {
    match panel {
        GeneticsPanel::Sequence => |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Sequence),
        GeneticsPanel::Alignment => |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Alignment),
        GeneticsPanel::Phylogenetics => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Phylogenetics)
        }
        GeneticsPanel::PopulationGenetics => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::PopulationGenetics)
        }
        GeneticsPanel::RnaStructure => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::RnaStructure)
        }
        GeneticsPanel::RnaDesigner => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::RnaDesigner)
        }
        GeneticsPanel::MolecularDynamics => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::MolecularDynamics)
        }
        GeneticsPanel::Cheminformatics => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Cheminformatics)
        }
        GeneticsPanel::MacromolecularStructure => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::MacromolecularStructure)
        }
        GeneticsPanel::QuantumChemistry => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::QuantumChemistry)
        }
        GeneticsPanel::Genomics => |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Genomics),
        GeneticsPanel::SystemsBiology => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::SystemsBiology)
        }
        GeneticsPanel::Docking => |app, ctx| draw_genetics(app, ctx, GeneticsPanel::Docking),
        GeneticsPanel::GeneEditing => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::GeneEditing)
        }
        GeneticsPanel::StructurePrediction => {
            |app, ctx| draw_genetics(app, ctx, GeneticsPanel::StructurePrediction)
        }
    }
}

fn draw_genetics(app: &mut ValenxApp, ctx: &egui::Context, panel: GeneticsPanel) {
    app.enable_genetics_workbench(panel);
    draw_genetics_workbench(app, ctx);
}

fn genetics_slug(panel: GeneticsPanel, index: usize) -> &'static str {
    // We need 'static slugs to plug into PanelSpec. Use a name table.
    match (panel, index) {
        (GeneticsPanel::Sequence, _) => "01-sequence",
        (GeneticsPanel::Alignment, _) => "02-alignment",
        (GeneticsPanel::Phylogenetics, _) => "03-phylogenetics",
        (GeneticsPanel::PopulationGenetics, _) => "04-population-genetics",
        (GeneticsPanel::RnaStructure, _) => "05-rna-structure",
        (GeneticsPanel::RnaDesigner, _) => "06-rna-designer",
        (GeneticsPanel::MolecularDynamics, _) => "07-molecular-dynamics",
        (GeneticsPanel::Cheminformatics, _) => "08-cheminformatics",
        (GeneticsPanel::MacromolecularStructure, _) => "09-macromolecular-structure",
        (GeneticsPanel::QuantumChemistry, _) => "10-quantum-chemistry",
        (GeneticsPanel::Genomics, _) => "11-genomics",
        (GeneticsPanel::SystemsBiology, _) => "12-systems-biology",
        (GeneticsPanel::Docking, _) => "13-docking",
        (GeneticsPanel::GeneEditing, _) => "14-gene-editing",
        (GeneticsPanel::StructurePrediction, _) => "15-structure-prediction",
    }
}

fn genetics_label(panel: GeneticsPanel) -> &'static str {
    match panel {
        GeneticsPanel::Sequence => "Sequence (bioseq)",
        GeneticsPanel::Alignment => "Alignment (align)",
        GeneticsPanel::Phylogenetics => "Phylogenetics (phylo)",
        GeneticsPanel::PopulationGenetics => "Population Genetics (popgen)",
        GeneticsPanel::RnaStructure => "RNA Structure (rnastruct)",
        GeneticsPanel::RnaDesigner => "RNA Designer (rnadesign)",
        GeneticsPanel::MolecularDynamics => "Molecular Dynamics (md)",
        GeneticsPanel::Cheminformatics => "Cheminformatics (cheminf)",
        GeneticsPanel::MacromolecularStructure => "Macromolecular Structure (biostruct)",
        GeneticsPanel::QuantumChemistry => "Quantum Chemistry (qchem)",
        GeneticsPanel::Genomics => "Genomics (genomics)",
        GeneticsPanel::SystemsBiology => "Systems Biology (sysbio)",
        GeneticsPanel::Docking => "Docking (dock-screen)",
        GeneticsPanel::GeneEditing => "Gene Editing (genediting)",
        GeneticsPanel::StructurePrediction => "Structure Prediction (structpredict)",
    }
}

fn genetics_description(panel: GeneticsPanel) -> &'static str {
    match panel {
        GeneticsPanel::Sequence => {
            "Sequence editing, translation, ORF finder, restriction \
             digest virtual gel, primer / PCR design (valenx-bioseq)."
        }
        GeneticsPanel::Alignment => {
            "Pairwise (Needleman-Wunsch / Smith-Waterman / Gotoh affine), \
             progressive MSA, k-mer seed search (valenx-align)."
        }
        GeneticsPanel::Phylogenetics => {
            "Distance-method phylogenetic trees (NJ / UPGMA), Newick \
             export, tree visualization (valenx-phylo)."
        }
        GeneticsPanel::PopulationGenetics => {
            "Wright-Fisher / coalescent simulations, summary stats (π, \
             θ, Tajima's D), allele-frequency trajectories (valenx-popgen)."
        }
        GeneticsPanel::RnaStructure => {
            "RNA secondary-structure folding (Zuker / LinearFold), \
             ensemble defect, MFE/PFE energy (valenx-rnastruct)."
        }
        GeneticsPanel::RnaDesigner => {
            "Guided synthetic-RNA design wizard — fold, visualize, \
             inverse design, mRNA LinearDesign, full construct \
             (valenx-rnadesign)."
        }
        GeneticsPanel::MolecularDynamics => {
            "Classical MD — Lennard-Jones, Coulomb, harmonic-bond + \
             angle + dihedral force fields, Berendsen / Langevin \
             integrators (valenx-md)."
        }
        GeneticsPanel::Cheminformatics => {
            "Cheminformatics descriptors (MW, logP, TPSA), Morgan \
             fingerprints, Tanimoto similarity, substructure search \
             (valenx-cheminf)."
        }
        GeneticsPanel::MacromolecularStructure => {
            "Macromolecular structure analysis — secondary structure, \
             Ramachandran plot, B-factor, contact map (valenx-biostruct)."
        }
        GeneticsPanel::QuantumChemistry => {
            "Hartree-Fock SCF, MP2 correlation energy, basis-set picker \
             (valenx-qchem)."
        }
        GeneticsPanel::Genomics => {
            "NGS / variant / CRISPR tooling — read alignment, variant \
             calling, gRNA design (valenx-genomics)."
        }
        GeneticsPanel::SystemsBiology => {
            "Systems-biology reaction networks — ODE / stochastic \
             simulation, sensitivity analysis (valenx-sysbio)."
        }
        GeneticsPanel::Docking => "Molecular docking + virtual screening (valenx-dock-screen).",
        GeneticsPanel::GeneEditing => {
            "CRISPR / base / prime / mRNA editing design (valenx-genediting)."
        }
        GeneticsPanel::StructurePrediction => {
            "Classical protein structure prediction (ab-initio fragment \
             assembly) + fixed-backbone design (valenx-structpredict)."
        }
    }
}

// ---------------------------------------------------------------------
// The off-screen wgpu + egui_wgpu render harness.
// ---------------------------------------------------------------------

struct ScreenshotHarness {
    device: wgpu::Device,
    queue: wgpu::Queue,
    target_format: wgpu::TextureFormat,
}

impl ScreenshotHarness {
    /// Build an off-screen wgpu device. Returns `None` if no adapter
    /// is available (CI sandbox without GPU + without a software
    /// fallback).
    fn new() -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
        }))?;
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("valenx.screenshots.device"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
            },
            None,
        ))
        .ok()?;
        eprintln!(
            "ScreenshotHarness: rendering on `{}` ({:?})",
            adapter.get_info().name,
            adapter.get_info().backend
        );
        // Linear (non-sRGB) target — the egui_wgpu renderer outputs
        // colour already gamma-encoded for its internal sRGB pipeline;
        // a non-sRGB attachment leaves those bytes alone so the PNG
        // matches what egui actually paints.
        Some(Self {
            device,
            queue,
            target_format: wgpu::TextureFormat::Rgba8Unorm,
        })
    }

    /// Render one panel into a PNG. Returns the texture dimensions so
    /// the caller can assert against them.
    fn render_panel_to_png(&mut self, spec: &PanelSpec, out: &Path) -> [u32; 2] {
        // --- Build the populated app + run one egui frame -----------
        let mut app = ValenxApp::default();
        let ctx = egui::Context::default();
        ctx.set_pixels_per_point(1.0);
        // Force every collapsing header / popup / menu open. Without
        // this, panels whose draw fn wraps their body in
        // `egui::CollapsingHeader::new(...).default_open(false)` (e.g.
        // every aero section) render only their one-line header bar
        // and the screenshot is useless for visual review.
        ctx.memory_mut(|m| m.set_everything_is_visible(true));
        // Bound the egui screen so SidePanels mount at a stable width.
        // Without this, the host SidePanel would extend off-screen.
        let raw = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::pos2(0.0, 0.0),
                egui::vec2(PANEL_WIDTH as f32, PANEL_HEIGHT as f32),
            )),
            ..Default::default()
        };
        let full_output = ctx.run(raw, |ctx| match spec.mount {
            Mount::HostSidePanel | Mount::CentralPanel => (spec.draw)(&mut app, ctx),
        });

        // --- Tessellate to triangle meshes ---------------------------
        let pixels_per_point = full_output.pixels_per_point;
        let paint_jobs = ctx.tessellate(full_output.shapes, pixels_per_point);
        let textures_delta = full_output.textures_delta;

        // --- Build the renderer + upload textures + buffers ---------
        let mut renderer = egui_wgpu::Renderer::new(
            &self.device,
            self.target_format,
            None, // no depth target
            1,    // no MSAA
        );
        for (id, image_delta) in &textures_delta.set {
            renderer.update_texture(&self.device, &self.queue, *id, image_delta);
        }
        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [PANEL_WIDTH, PANEL_HEIGHT],
            pixels_per_point,
        };

        // --- Off-screen color target + readback buffer --------------
        let color_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("valenx.screenshots.color"),
            size: wgpu::Extent3d {
                width: PANEL_WIDTH,
                height: PANEL_HEIGHT,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: self.target_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // 256-byte-aligned row pitch (wgpu COPY_BYTES_PER_ROW_ALIGNMENT).
        let bytes_per_pixel = 4u32;
        let unpadded_bytes_per_row = PANEL_WIDTH * bytes_per_pixel;
        let padded_bytes_per_row = align_up(unpadded_bytes_per_row, 256);
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("valenx.screenshots.readback"),
            size: (padded_bytes_per_row * PANEL_HEIGHT) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        // --- Encode + submit ----------------------------------------
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("valenx.screenshots.encoder"),
            });
        let user_cmd_buffers = renderer.update_buffers(
            &self.device,
            &self.queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );
        {
            let mut rpass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("valenx.screenshots.pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        // egui's panel-frame fill will paint the
                        // background on its own; a clear colour here
                        // shouldn't be visible. Pick the egui default
                        // dark panel background as a sanity fallback
                        // (RGB ~ 27, 27, 27 / 255).
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: 0.106,
                            g: 0.106,
                            b: 0.106,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            renderer.render(&mut rpass, &paint_jobs, &screen_descriptor);
        }
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &readback,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(padded_bytes_per_row),
                    rows_per_image: Some(PANEL_HEIGHT),
                },
            },
            wgpu::Extent3d {
                width: PANEL_WIDTH,
                height: PANEL_HEIGHT,
                depth_or_array_layers: 1,
            },
        );

        // Submit the user command buffers (texture upload prep) FIRST
        // — they must happen before the render pass / copy reads from
        // those buffers.
        let mut cmd_buffers = user_cmd_buffers;
        cmd_buffers.push(encoder.finish());
        self.queue.submit(cmd_buffers);

        // --- Free egui textures we won't reuse ----------------------
        for id in &textures_delta.free {
            renderer.free_texture(id);
        }

        // --- Map + read back pixels ---------------------------------
        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .expect("map_async callback never fired")
            .expect("readback map failed");
        let data = slice.get_mapped_range();

        // Strip the row padding back to a tight unpadded image.
        let mut tight: Vec<u8> =
            Vec::with_capacity((unpadded_bytes_per_row * PANEL_HEIGHT) as usize);
        for row in 0..PANEL_HEIGHT {
            let start = (row * padded_bytes_per_row) as usize;
            let end = start + unpadded_bytes_per_row as usize;
            tight.extend_from_slice(&data[start..end]);
        }
        drop(data);
        readback.unmap();

        // --- Assert the frame is not uniformly the clear colour -----
        // The clear colour above was (27, 27, 27, 255) — any frame that
        // came out 100 % identical to that means the render produced
        // nothing. The threshold is deliberately conservative (~100 px,
        // < 0.01 % of frame) — even a collapsed CollapsingHeader paints
        // hundreds of glyph pixels for its label, so anything below
        // this floor means the draw fn truly no-op'd.
        let clear_rgb = [27u8, 27u8, 27u8];
        let non_clear = tight
            .chunks_exact(4)
            .filter(|p| [p[0], p[1], p[2]] != clear_rgb)
            .count();
        let total_px = (PANEL_WIDTH * PANEL_HEIGHT) as usize;
        assert!(
            non_clear > 100,
            "panel `{}` rendered as a near-uniform clear-colour frame \
             — only {non_clear}/{total_px} pixels differ from the clear \
             colour; the draw fn likely no-op'd",
            spec.slug,
        );

        // --- Encode + write PNG -------------------------------------
        write_png(out, &tight, PANEL_WIDTH, PANEL_HEIGHT)
            .unwrap_or_else(|e| panic!("write_png({}): {e}", out.display()));

        [PANEL_WIDTH, PANEL_HEIGHT]
    }
}

fn align_up(value: u32, alignment: u32) -> u32 {
    (value + alignment - 1) & !(alignment - 1)
}

fn write_png(out: &Path, rgba: &[u8], width: u32, height: u32) -> Result<(), String> {
    let file = fs::File::create(out).map_err(|e| e.to_string())?;
    let w = std::io::BufWriter::new(file);
    let mut encoder = png::Encoder::new(w, width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().map_err(|e| e.to_string())?;
    writer.write_image_data(rgba).map_err(|e| e.to_string())?;
    Ok(())
}

// ---------------------------------------------------------------------
// Workspace-root resolution + INDEX.md emission
// ---------------------------------------------------------------------

/// Resolve the workspace root from `CARGO_MANIFEST_DIR` (which points
/// at `crates/valenx-app` when running `cargo test -p valenx-app`).
fn workspace_root() -> PathBuf {
    // CARGO_MANIFEST_DIR -> crates/valenx-app; workspace root is two
    // levels up.
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest_dir
        .parent()
        .and_then(Path::parent)
        .unwrap_or(&manifest_dir)
        .to_path_buf()
}

/// Write a `screenshots/INDEX.md` summarising every PNG.
fn write_index(out_root: &Path, rendered: &[RenderedPanel]) -> std::io::Result<()> {
    let mut by_workbench: BTreeMap<&'static str, Vec<&RenderedPanel>> = BTreeMap::new();
    for r in rendered {
        by_workbench
            .entry(r.spec.workbench.label())
            .or_default()
            .push(r);
    }

    let mut md = String::new();
    md.push_str("# Valenx workbench screenshots\n\n");
    md.push_str(
        "Generated by `cargo test -p valenx-app --test \
         headless_screenshots`. One PNG per workbench panel, captured \
         in a representative populated state.\n\n",
    );
    md.push_str(&format!(
        "Image size: **{PANEL_WIDTH} × {PANEL_HEIGHT}** \
         (`Rgba8Unorm`), one frame per panel.\n\n",
    ));

    for (workbench, panels) in &by_workbench {
        md.push_str(&format!("## {workbench}\n\n"));
        md.push_str("| Panel | Description | PNG |\n");
        md.push_str("|---|---|---|\n");
        for p in panels {
            let rel = p
                .out_path
                .strip_prefix(out_root)
                .unwrap_or(&p.out_path)
                .display()
                .to_string()
                .replace('\\', "/");
            md.push_str(&format!(
                "| **{}** | {} | [`{}`]({}) |\n",
                p.spec.label, p.spec.description, rel, rel,
            ));
        }
        md.push('\n');
    }

    fs::write(out_root.join("INDEX.md"), md)
}
