//! "File → New Project…" modal dialog.
//!
//! Surfaces the same scaffolding API the `valenx-init` CLI exposes —
//! [`valenx_core::init_templates::scaffold_project`] — so the GUI no
//! longer requires the user to drop down to a terminal just to create
//! a fresh project on first launch.
//!
//! The dialog is intentionally read-only with respect to disk until
//! the user clicks **Create**. The folder picker (rfd) is the only
//! native dialog reached from here, and it's gated behind an explicit
//! button click so the test harness never blocks on it.
//!
//! ## Flow
//!
//! 1. User opens the dialog (File menu, command palette, or Ctrl+N).
//! 2. Picks a name, a template, and a target directory. The default
//!    directory is `~/Documents/Valenx Projects/<name>/` so the
//!    project lands somewhere users routinely browse.
//! 3. Clicks **Create**. We validate the inputs, call
//!    [`scaffold_project`], and on success trigger the existing
//!    `ValenxApp::load_project` flow so the project appears in the
//!    Browser tree.
//! 4. On any validation / scaffold error the dialog stays open with
//!    the error rendered inline.

use std::path::PathBuf;

use eframe::egui;

use valenx_core::init_templates::{scaffold_project, Template};

/// One entry in the template chooser. Pairs a [`Template`] with the
/// display name + one-line description the chooser renders. The
/// chooser groups these by [`TemplateGroup`] so the dropdown isn't a
/// 150-item flat list.
#[derive(Clone, Copy, Debug)]
pub struct TemplateOption {
    pub template: Template,
    pub display_name: &'static str,
    pub description: &'static str,
}

/// Top-level grouping for the chooser. Each group renders as a
/// collapsing header in the dialog. The order here is the order the
/// groups appear; the most common categories (Empty, Engineering,
/// Physics & Chemistry) come first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TemplateGroup {
    Empty,
    Engineering,
    PhysicsChemistry,
    Bioinformatics,
    StructurePrediction,
    Workflows,
}

impl TemplateGroup {
    /// Display label used as the collapsing-header text.
    pub fn label(self) -> &'static str {
        match self {
            Self::Empty => "Start from scratch",
            Self::Engineering => "Engineering simulations",
            Self::PhysicsChemistry => "Physics & chemistry",
            Self::Bioinformatics => "Bioinformatics",
            Self::StructurePrediction => "Protein structure prediction",
            Self::Workflows => "Workflow managers",
        }
    }

    /// Iteration order for the chooser (most common first).
    pub const ALL: [Self; 6] = [
        Self::Empty,
        Self::Engineering,
        Self::PhysicsChemistry,
        Self::Bioinformatics,
        Self::StructurePrediction,
        Self::Workflows,
    ];
}

/// The full curated chooser catalogue. Each row carries its parent
/// group so the dialog can render group-by-group without a parallel
/// lookup table.
///
/// Kept much smaller than `init_templates::template_rows()` (which
/// has 150+ entries) — the chooser surfaces a representative subset
/// covering the major adapters in each domain. Users who need a more
/// obscure template can still drop down to `valenx-init`.
pub fn template_catalogue() -> &'static [(TemplateGroup, TemplateOption)] {
    CATALOGUE
}

#[rustfmt::skip]
const CATALOGUE: &[(TemplateGroup, TemplateOption)] = &[
    // -- Empty --
    (TemplateGroup::Empty, TemplateOption {
        template: Template::Empty,
        display_name: "Empty project",
        description: "Minimal skeleton — no per-physics block. Edit case.toml to bring your own solver.",
    }),

    // -- Engineering simulations --
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::Cfd,
        display_name: "CFD (OpenFOAM cavity)",
        description: "OpenFOAM simpleFoam — incompressible RANS over a lid-driven cavity.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::Su2,
        display_name: "SU2 NACA0012",
        description: "SU2 compressible CFD — NACA 0012 airfoil starter.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::Fea,
        display_name: "FEA (CalculiX cantilever)",
        description: "CalculiX linear-static — tip-loaded cantilever beam.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::CodeAster,
        display_name: "Static beam (Code_Aster)",
        description: "Code_Aster `as_run` on a user-built `.export`.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::OpenRadioss,
        display_name: "Drop test (OpenRadioss)",
        description: "OpenRadioss explicit dynamics — engine-deck only.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::Netgen,
        display_name: "CSG box (Netgen)",
        description: "Netgen CSG meshing — axis-aligned unit cube.",
    }),
    (TemplateGroup::Engineering, TemplateOption {
        template: Template::Gmsh,
        display_name: "Box mesh (gmsh)",
        description: "gmsh procedural meshing — Delaunay tet box.",
    }),

    // -- Physics & chemistry --
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Chemistry,
        display_name: "Chemistry (Cantera CH4/air)",
        description: "Cantera equilibrium-HP — methane / air mixture.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::ElmerHeat,
        display_name: "Heat conduction (Elmer)",
        description: "Elmer steady heat — two pinned-temperature faces.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Meep,
        display_name: "Ring resonator (Meep)",
        description: "Meep FDTD — Python ring-resonator script.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Lammps,
        display_name: "Lennard-Jones (LAMMPS)",
        description: "LAMMPS classical MD — Lennard-Jones FCC fluid (NVE).",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Gromacs,
        display_name: "Lysozyme MD (GROMACS)",
        description: "GROMACS `gmx mdrun` on a user-built `.tpr`.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Openmm,
        display_name: "Protein relax (OpenMM)",
        description: "OpenMM — Python-native MD minimisation + DCD output.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Psi4,
        display_name: "Quantum chem (Psi4)",
        description: "Psi4 — HF/DFT/post-HF quantum chemistry.",
    }),
    (TemplateGroup::PhysicsChemistry, TemplateOption {
        template: Template::Xtb,
        display_name: "Semiempirical QC (xTB)",
        description: "xTB — extended tight-binding semiempirical quantum chemistry.",
    }),

    // -- Bioinformatics --
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Biopython,
        display_name: "Biopython analyse",
        description: "Biopython — sequence / structural-bio Python script.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Rdkit,
        display_name: "RDKit screen",
        description: "RDKit — cheminformatics Python script.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Chimerax,
        display_name: "ChimeraX render",
        description: "ChimeraX — `.cxc` command-script renderer.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Oxdna,
        display_name: "oxDNA duplex",
        description: "oxDNA — coarse-grained DNA / RNA MD on `input.dat`.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Mdanalysis,
        display_name: "MDAnalysis trajectory RMSD",
        description: "MDAnalysis — trajectory analysis Python script.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Bwa,
        display_name: "BWA align (short-read)",
        description: "BWA — short-read DNA alignment via `bwa mem`.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Minimap2,
        display_name: "Minimap2 align (long-read)",
        description: "minimap2 — long-read + cross-domain alignment.",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Samtools,
        display_name: "samtools flagstat",
        description: "samtools — SAM/BAM multitool (view / sort / index / flagstat).",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Bcftools,
        display_name: "Variant calling (bcftools)",
        description: "bcftools — VCF/BCF multitool (view / call / filter / concat).",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::ViennaRna,
        display_name: "ViennaRNA fold",
        description: "ViennaRNA RNAfold — secondary-structure prediction (academic).",
    }),
    (TemplateGroup::Bioinformatics, TemplateOption {
        template: Template::Scanpy,
        display_name: "Scanpy single-cell",
        description: "Scanpy — Python single-cell analysis.",
    }),

    // -- Structure prediction --
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::Colabfold,
        display_name: "ColabFold folding",
        description: "ColabFold — protein structure prediction from FASTA.",
    }),
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::Esmfold,
        display_name: "ESMFold predict",
        description: "ESMFold — Meta protein language model structure prediction.",
    }),
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::Alphafold2,
        display_name: "AlphaFold 2 predict",
        description: "AlphaFold 2 — DeepMind structure prediction (open weights).",
    }),
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::Alphafold3,
        display_name: "AlphaFold 3 predict",
        description: "AlphaFold 3 — all-atom complex prediction (non-commercial weights).",
    }),
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::OmegaFold,
        display_name: "OmegaFold predict",
        description: "OmegaFold — single-sequence structure prediction (no MSA).",
    }),
    (TemplateGroup::StructurePrediction, TemplateOption {
        template: Template::Foldseek,
        display_name: "FoldSeek search",
        description: "FoldSeek — protein structure search via 3Di alphabet.",
    }),

    // -- Workflow managers --
    (TemplateGroup::Workflows, TemplateOption {
        template: Template::Nextflow,
        display_name: "Nextflow pipeline",
        description: "Nextflow — pipeline orchestrator.",
    }),
    (TemplateGroup::Workflows, TemplateOption {
        template: Template::Snakemake,
        display_name: "Snakemake pipeline",
        description: "Snakemake — rule-based pipeline orchestrator.",
    }),
];

/// In-flight modal state. `None` on `ValenxApp` means the dialog is
/// closed. The render call mutates `name` / `location` / `template` /
/// `error` in-place every frame.
#[derive(Clone, Debug)]
pub struct NewProjectDialog {
    /// Project name — populates `[project] name = "..."` in the
    /// generated `project.toml` and (when the location field is in
    /// auto-track mode) drives the default folder name too.
    pub name: String,
    /// Target directory. Defaults to
    /// `~/Documents/Valenx Projects/<name>/`; the user can override
    /// via the Browse button.
    pub location: String,
    /// When `true`, edits to `name` re-derive `location`. Flipped to
    /// `false` the first time the user picks a folder explicitly so
    /// their pick isn't clobbered by later name changes.
    pub auto_location: bool,
    /// Currently selected template.
    pub template: Template,
    /// User-visible error after a failed validation / scaffold. Cleared
    /// on the next Create attempt.
    pub error: Option<String>,
}

impl Default for NewProjectDialog {
    fn default() -> Self {
        let name = String::from("my-project");
        let location = default_project_location(&name);
        Self {
            name,
            location,
            auto_location: true,
            template: Template::Empty,
            error: None,
        }
    }
}

/// Derive the default project location for a given project name.
/// Returns `~/Documents/Valenx Projects/<name>/` when HOME-equivalent
/// resolves, otherwise the system temp dir as a sensible fallback so
/// the field never starts empty.
pub fn default_project_location(name: &str) -> String {
    let base = documents_dir().unwrap_or_else(std::env::temp_dir);
    base.join("Valenx Projects")
        .join(name)
        .to_string_lossy()
        .to_string()
}

/// Resolve the per-OS Documents folder. Picks `~/Documents` on
/// Windows + macOS + Linux (the freedesktop standard); returns `None`
/// if no HOME-equivalent env var is set.
fn documents_dir() -> Option<PathBuf> {
    if cfg!(target_os = "windows") {
        std::env::var_os("USERPROFILE")
            .map(|p| PathBuf::from(p).join("Documents"))
            .or_else(|| std::env::var_os("HOME").map(|p| PathBuf::from(p).join("Documents")))
    } else {
        std::env::var_os("HOME").map(|p| PathBuf::from(p).join("Documents"))
    }
}

/// Validate a project name. Rejects empty, whitespace-only, and any
/// name containing characters that aren't legal on every supported
/// filesystem (Windows-9x reserved punctuation is the strictest, so
/// we use that set everywhere for portability).
pub fn validate_name(name: &str) -> Result<(), String> {
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return Err("Project name can't be empty.".into());
    }
    // Windows-reserved chars plus the slash variants that confuse
    // every loader. Mirror cmd's "you can't name a file this" list so
    // the project lands on every supported OS.
    const BAD: &[char] = &['/', '\\', ':', '*', '?', '"', '<', '>', '|'];
    if let Some(bad) = trimmed.chars().find(|c| BAD.contains(c)) {
        return Err(format!(
            "Project name can't contain `{bad}` (illegal on Windows / Linux / macOS)."
        ));
    }
    Ok(())
}

/// What [`render`] tells the caller to do once it returns. The
/// dialog has no app-mutating side effects of its own — every disk
/// action funnels through this enum so callers control the post-
/// create project-load flow.
#[derive(Clone, Debug)]
pub enum DialogOutcome {
    /// Nothing changed this frame; keep the dialog open.
    Stay,
    /// User clicked Cancel (or hit Esc). Close the dialog, do not
    /// load anything.
    Cancel,
    /// User clicked Create and validation + scaffolding succeeded.
    /// Caller should close the dialog and call
    /// `app.load_project(<path>)`.
    Created(PathBuf),
}

/// Render the New Project modal. Returns the outcome the caller acts
/// on next.
///
/// The function mutates `dialog` in-place every frame (text input
/// state, template pick, error display). It only triggers a disk
/// write when the user clicks **Create** — folder picker dialogs are
/// gated behind explicit button clicks.
pub fn render(ctx: &egui::Context, dialog: &mut NewProjectDialog) -> DialogOutcome {
    let mut outcome = DialogOutcome::Stay;
    let mut keep_open = true;
    let mut create_clicked = false;
    let mut cancel_clicked = false;
    let mut browse_clicked = false;
    let mut name_changed = false;

    egui::Window::new("New Project")
        .open(&mut keep_open)
        .collapsible(false)
        .resizable(true)
        .default_size([640.0, 540.0])
        .min_width(520.0)
        .show(ctx, |ui| {
            ui.add_space(4.0);
            ui.heading("Create a new Valenx project");
            ui.label(
                "Pick a name, a template, and where the new `.valenx` directory \
                 should land. Valenx will scaffold the project.toml + a starter \
                 case so you can run something straight away.",
            );
            ui.add_space(10.0);

            // Project name field — labelled to align with the next two rows.
            ui.horizontal(|ui| {
                ui.label("Name:");
                let edit = egui::TextEdit::singleline(&mut dialog.name)
                    .hint_text("my-project")
                    .desired_width(360.0);
                if ui.add(edit).changed() {
                    name_changed = true;
                }
            });

            // Location row — text field + Browse button. The
            // text field is editable so users can paste a path
            // without going through the picker.
            ui.horizontal(|ui| {
                ui.label("Location:");
                let edit = egui::TextEdit::singleline(&mut dialog.location)
                    .desired_width(360.0)
                    .hint_text("(picks Documents/Valenx Projects/<name> by default)");
                if ui.add(edit).changed() {
                    // Manual edit breaks auto-track so the user's
                    // path isn't clobbered by subsequent name edits.
                    dialog.auto_location = false;
                }
                if ui.button("Browse…").clicked() {
                    browse_clicked = true;
                }
            });
            ui.add_space(8.0);

            // Template chooser — grouped collapsing headers each
            // containing a selectable list. The "Empty" group opens
            // by default so the first-launch user sees a usable
            // chooser without scrolling.
            ui.label(egui::RichText::new("Template").strong());
            ui.add_space(2.0);
            egui::ScrollArea::vertical()
                .max_height(280.0)
                .show(ui, |ui| {
                    for group in TemplateGroup::ALL {
                        let default_open = matches!(
                            group,
                            TemplateGroup::Empty | TemplateGroup::Engineering
                        );
                        egui::CollapsingHeader::new(group.label())
                            .default_open(default_open)
                            .show(ui, |ui| {
                                for (g, opt) in template_catalogue() {
                                    if *g != group {
                                        continue;
                                    }
                                    // `Template` derives `PartialEq + Eq`, so
                                    // `==` is the direct way to express this.
                                    // The previous `discriminant`-based check
                                    // was a workaround from before the derive
                                    // landed; it accidentally treated payload-
                                    // bearing variants as equal across their
                                    // payloads. Direct equality is correct
                                    // for the current Template enum which
                                    // contains only unit variants today.
                                    let selected = dialog.template == opt.template;
                                    let label = format!(
                                        "{}  —  {}",
                                        opt.display_name, opt.description
                                    );
                                    if ui
                                        .selectable_label(selected, label)
                                        .on_hover_text(opt.description)
                                        .clicked()
                                    {
                                        dialog.template = opt.template;
                                    }
                                }
                            });
                    }
                });

            // Error display — only painted when set, so the dialog
            // height doesn't oscillate between renders. Pinned just
            // above the action row.
            if let Some(err) = &dialog.error {
                ui.add_space(8.0);
                ui.colored_label(egui::Color32::from_rgb(220, 80, 80), err);
            }

            ui.add_space(10.0);
            ui.horizontal(|ui| {
                if ui
                    .add(egui::Button::new(
                        egui::RichText::new("Create").strong(),
                    ))
                    .clicked()
                {
                    create_clicked = true;
                }
                if ui.button("Cancel").clicked() {
                    cancel_clicked = true;
                }
            });
        });

    // Window close button (the [x] on the title bar) flips
    // `keep_open` to false. Treat that the same as Cancel.
    if !keep_open {
        cancel_clicked = true;
    }

    if cancel_clicked {
        return DialogOutcome::Cancel;
    }

    // Auto-track location to name edits as long as the user hasn't
    // manually picked a path. Without this the default would freeze
    // at "my-project" even when the user types a real name.
    if name_changed && dialog.auto_location {
        dialog.location = default_project_location(dialog.name.trim());
    }

    if browse_clicked {
        // The folder picker is a native modal that blocks the main
        // thread until the user dismisses it. We seed it with the
        // current location (or its parent) so the picker opens at
        // a meaningful starting point.
        let seed = PathBuf::from(&dialog.location);
        let start = if seed.exists() {
            seed
        } else {
            seed.parent().map(PathBuf::from).unwrap_or_default()
        };
        let mut picker = rfd::FileDialog::new().set_title("Choose project location");
        if start.exists() {
            picker = picker.set_directory(start);
        }
        if let Some(picked) = picker.pick_folder() {
            dialog.location = picked.to_string_lossy().to_string();
            dialog.auto_location = false;
        }
    }

    if create_clicked {
        match try_create(dialog) {
            Ok(path) => {
                outcome = DialogOutcome::Created(path);
            }
            Err(reason) => {
                dialog.error = Some(reason);
            }
        }
    }

    outcome
}

/// Validate inputs and call into [`scaffold_project`]. Pure
/// function: takes `dialog` immutably + does no UI.
///
/// On success returns the resolved project path the caller should
/// load. On failure returns a user-facing reason the dialog renders
/// inline.
pub fn try_create(dialog: &NewProjectDialog) -> Result<PathBuf, String> {
    validate_name(&dialog.name)?;
    let trimmed_name = dialog.name.trim().to_string();
    let location = dialog.location.trim();
    if location.is_empty() {
        return Err("Location can't be empty — pick a folder.".into());
    }
    let target = PathBuf::from(location);
    // Refuse to scaffold into a directory that already contains a
    // populated project — mirrors the CLI's behaviour so users
    // don't accidentally clobber existing work.
    if target.join("project.toml").exists() {
        return Err(format!(
            "{} already contains a project.toml — pick an empty location or delete it first.",
            target.display()
        ));
    }
    // The parent dir is what scaffold_project will create under, so
    // we make sure it's a writable real directory before we attempt
    // the scaffold. The leaf (target itself) is created by
    // scaffold_project via create_dir_all so we don't need to
    // pre-create it.
    if let Some(parent) = target.parent() {
        if !parent.exists() {
            // Try to create the parent — this covers the common case
            // of the default `~/Documents/Valenx Projects/` not
            // existing yet on a fresh OS install.
            std::fs::create_dir_all(parent).map_err(|e| {
                format!(
                    "Couldn't create parent directory {}: {e}",
                    parent.display()
                )
            })?;
        } else if !parent.is_dir() {
            return Err(format!(
                "Parent path {} exists but isn't a directory.",
                parent.display()
            ));
        }
    }
    scaffold_project(&target, dialog.template, Some(&trimmed_name))
        .map_err(|e| format!("Couldn't create project: {e}"))?;
    Ok(target)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn validate_name_rejects_empty_and_whitespace() {
        assert!(validate_name("").is_err());
        assert!(validate_name("   ").is_err());
        assert!(validate_name("\t\n").is_err());
    }

    #[test]
    fn validate_name_rejects_path_separators_and_reserved_chars() {
        // Windows-illegal characters that should be rejected on every
        // OS so the project lands cleanly when copied across systems.
        for bad in ['/', '\\', ':', '*', '?', '"', '<', '>', '|'] {
            let name = format!("bad{bad}name");
            assert!(
                validate_name(&name).is_err(),
                "expected {name:?} to be rejected"
            );
        }
    }

    #[test]
    fn validate_name_accepts_alphanumeric_and_hyphen() {
        assert!(validate_name("my-project").is_ok());
        assert!(validate_name("project_42").is_ok());
        assert!(validate_name("CFDStudy.v2").is_ok());
    }

    #[test]
    fn default_dialog_has_sensible_defaults() {
        let d = NewProjectDialog::default();
        assert_eq!(d.name, "my-project");
        assert!(matches!(d.template, Template::Empty));
        assert!(d.auto_location);
        assert!(!d.location.is_empty());
        assert!(d.error.is_none());
    }

    #[test]
    fn template_catalogue_is_non_empty_and_covers_every_group() {
        let cat = template_catalogue();
        assert!(!cat.is_empty());
        let mut groups_seen = std::collections::HashSet::new();
        for (g, _) in cat {
            groups_seen.insert(*g);
        }
        for g in TemplateGroup::ALL {
            assert!(
                groups_seen.contains(&g),
                "no template surfaces under group {g:?}"
            );
        }
    }

    #[test]
    fn try_create_scaffolds_into_a_clean_dir() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-new-project-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        // Don't pre-create — scaffold_project handles it.
        let dialog = NewProjectDialog {
            name: "smoke-test".into(),
            location: tmp.to_string_lossy().to_string(),
            auto_location: false,
            template: Template::Empty,
            error: None,
        };
        let result = try_create(&dialog);
        assert!(result.is_ok(), "expected success, got {result:?}");
        let path = result.unwrap();
        assert_eq!(path, tmp);
        assert!(tmp.join("project.toml").is_file(), "project.toml missing");
        assert!(tmp.join("cases").is_dir(), "cases dir missing");
        // Cleanup — best-effort, don't fail the test on cleanup errors.
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn try_create_refuses_to_overwrite_populated_dir() {
        let tmp = std::env::temp_dir().join(format!(
            "valenx-new-project-overwrite-test-{}",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        std::fs::create_dir_all(&tmp).expect("tmp dir");
        std::fs::write(tmp.join("project.toml"), "stub").expect("seed");
        let dialog = NewProjectDialog {
            name: "smoke-test".into(),
            location: tmp.to_string_lossy().to_string(),
            auto_location: false,
            template: Template::Empty,
            error: None,
        };
        let result = try_create(&dialog);
        assert!(result.is_err(), "expected refusal, got {result:?}");
        let _ = std::fs::remove_dir_all(&tmp);
    }

    #[test]
    fn try_create_rejects_invalid_name() {
        let dialog = NewProjectDialog {
            name: "bad/name".into(),
            location: std::env::temp_dir().to_string_lossy().to_string(),
            auto_location: false,
            template: Template::Empty,
            error: None,
        };
        let result = try_create(&dialog);
        assert!(result.is_err(), "expected rejection, got {result:?}");
    }

    #[test]
    fn default_project_location_includes_name() {
        let loc = default_project_location("foo");
        assert!(loc.ends_with("foo") || loc.ends_with("foo/") || loc.ends_with("foo\\"));
        assert!(loc.contains("Valenx Projects"));
    }
}
