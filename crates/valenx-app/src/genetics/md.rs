//! Panel 6 — **Molecular Dynamics** (`valenx-md`).
//!
//! Load a structure (PDB or XYZ), build a real molecular-mechanics
//! force field for it, configure and run a classical molecular-
//! dynamics simulation, and inspect the energy / temperature traces —
//! all native `valenx-md` calls.
//!
//! `valenx-md`'s PDB / XYZ readers produce a *bond-free* topology
//! (PDB `CONECT` is not parsed, XYZ has no connectivity at all). This
//! panel closes that gap: it detects bonds from interatomic distances
//! with the covalent-radius rule, derives the full bonded topology
//! (bonds → angles → proper dihedrals), and parameterises it with a
//! small generic harmonic force field plus **element-specific**
//! Lennard-Jones parameters (a built-in table rather than one global
//! σ/ε). The user picks the integrator, the thermostat, the step
//! count and the temperature; initial velocities are drawn from the
//! Maxwell-Boltzmann distribution at that temperature. After the run
//! the panel plots kinetic / potential / total energy and temperature
//! and reports the Cα-free all-atom RMSD from the starting structure.
//!
//! The built structure can be pushed into the app's 3-D viewport as a
//! ball-and-stick model via [`crate::genetics::molecule_view`].

use eframe::egui;
use egui_plot::{Legend, Line, PlotPoints};
use nalgebra::Vector3;

use valenx_md::analysis::rmsd::rmsd as kabsch_rmsd;
use valenx_md::ensemble::andersen::{Andersen, VelocityRescale};
use valenx_md::ensemble::berendsen::Berendsen;
use valenx_md::ensemble::Thermostat;
use valenx_md::forcefield::{
    AngleParam, BondParam, CombiningRule, DihedralParam, ForceField, LjParam,
};
use valenx_md::integrate::leapfrog::LeapFrog;
use valenx_md::integrate::velocity_verlet::VelocityVerlet;
use valenx_md::integrate::Integrator;
use valenx_md::io::pdb::read_pdb;
use valenx_md::io::xyz::read_xyz;
use valenx_md::rng::Rng;
use valenx_md::sim::Simulation;
use valenx_md::system::System;

use super::common;
use super::molecule_view::{self, ViewMolecule};
use crate::plot_ui::managed_plot_mem_cfg;
use crate::ValenxApp;

/// Which time integrator the run uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum IntegratorChoice {
    /// Velocity-Verlet — the MD workhorse, exactly time-reversible.
    #[default]
    VelocityVerlet,
    /// Leapfrog — the GROMACS default, half-step velocities.
    Leapfrog,
}

impl IntegratorChoice {
    fn label(self) -> &'static str {
        match self {
            IntegratorChoice::VelocityVerlet => "Velocity-Verlet",
            IntegratorChoice::Leapfrog => "Leapfrog",
        }
    }
}

/// Which temperature-coupling scheme (thermostat) the run uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum ThermostatChoice {
    /// NVE — no thermostat, energy-conserving.
    #[default]
    None,
    /// Berendsen weak coupling — fast equilibration, not a true
    /// canonical ensemble.
    Berendsen,
    /// Andersen — stochastic collisions; a true canonical ensemble.
    Andersen,
    /// Velocity-rescale (Bussi) — canonical ensemble, smooth.
    VelocityRescale,
}

impl ThermostatChoice {
    fn label(self) -> &'static str {
        match self {
            ThermostatChoice::None => "None (NVE)",
            ThermostatChoice::Berendsen => "Berendsen",
            ThermostatChoice::Andersen => "Andersen",
            ThermostatChoice::VelocityRescale => "Velocity-rescale (Bussi)",
        }
    }
}

/// Snapshot of every editable input the MD panel owns. Numerics +
/// the structure text + the engine choices — `Ctrl+Z` rewinds them
/// all atomically.
#[derive(Clone, Debug, PartialEq, Default)]
pub(crate) struct MdSnapshot {
    pub(crate) structure_text: String,
    pub(crate) is_xyz: bool,
    pub(crate) bond_tolerance: f64,
    pub(crate) integrator: IntegratorChoice,
    pub(crate) thermostat: ThermostatChoice,
    pub(crate) temperature: f64,
    pub(crate) timestep_fs: f64,
    pub(crate) steps: usize,
    pub(crate) report_interval: usize,
    pub(crate) seed: u64,
}

/// Form + result state for the Molecular Dynamics panel.
pub struct MdPanel {
    /// Raw structure text (PDB or XYZ).
    structure_text: String,
    /// `true` when the loaded text is XYZ, `false` for PDB.
    is_xyz: bool,
    /// Multiplier on the covalent-radius bond-detection cutoff. `1.0`
    /// is the textbook rule; a small slack absorbs coordinate noise.
    bond_tolerance: f64,
    /// Integrator selection.
    integrator: IntegratorChoice,
    /// Thermostat selection.
    thermostat: ThermostatChoice,
    /// Target / initial temperature (K).
    temperature: f64,
    /// Integration time step (fs).
    timestep_fs: f64,
    /// Number of integration steps.
    steps: usize,
    /// State-report interval.
    report_interval: usize,
    /// RNG seed for the Maxwell-Boltzmann velocity draw + stochastic
    /// thermostats.
    seed: u64,
    error: Option<String>,
    result: String,
    /// Energy trace: `(step, potential, kinetic, total)`.
    trace: Vec<(f64, f64, f64, f64)>,
    /// Temperature trace: `(step, temperature)`.
    temp_trace: Vec<(f64, f64)>,
    /// Undo / redo for every input the panel exposes.
    history: crate::undo::History<MdSnapshot>,
}

impl MdPanel {
    fn snapshot(&self) -> MdSnapshot {
        MdSnapshot {
            structure_text: self.structure_text.clone(),
            is_xyz: self.is_xyz,
            bond_tolerance: self.bond_tolerance,
            integrator: self.integrator,
            thermostat: self.thermostat,
            temperature: self.temperature,
            timestep_fs: self.timestep_fs,
            steps: self.steps,
            report_interval: self.report_interval,
            seed: self.seed,
        }
    }
    fn restore(&mut self, s: MdSnapshot) {
        self.structure_text = s.structure_text;
        self.is_xyz = s.is_xyz;
        self.bond_tolerance = s.bond_tolerance;
        self.integrator = s.integrator;
        self.thermostat = s.thermostat;
        self.temperature = s.temperature;
        self.timestep_fs = s.timestep_fs;
        self.steps = s.steps;
        self.report_interval = s.report_interval;
        self.seed = s.seed;
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

impl Default for MdPanel {
    fn default() -> Self {
        MdPanel {
            // A small ethane-like XYZ — eight atoms with a C-C bond,
            // six C-H bonds, real angles and dihedrals once the bond
            // detector runs. Folds into a genuine bonded MD run with
            // no file I/O.
            structure_text: "8\nethane\n\
                             C  -0.762  0.000  0.000\n\
                             C   0.762  0.000  0.000\n\
                             H  -1.156  0.512  0.886\n\
                             H  -1.156  0.512 -0.886\n\
                             H  -1.156 -1.025  0.000\n\
                             H   1.156 -0.512  0.886\n\
                             H   1.156 -0.512 -0.886\n\
                             H   1.156  1.025  0.000\n"
                .to_string(),
            is_xyz: true,
            bond_tolerance: 1.0,
            integrator: IntegratorChoice::VelocityVerlet,
            thermostat: ThermostatChoice::Berendsen,
            temperature: 300.0,
            timestep_fs: 1.0,
            steps: 400,
            report_interval: 10,
            seed: 1,
            error: None,
            result: String::new(),
            trace: Vec::new(),
            temp_trace: Vec::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Element-specific Lennard-Jones parameters, GROMACS units (σ in nm,
/// ε in kJ/mol). A small AMBER/OPLS-flavoured table covering the
/// elements common in biomolecular structures; an unknown type falls
/// back to a carbon-like atom so a force-field build never fails.
fn element_lj(element: &str) -> LjParam {
    let (sigma, epsilon) = match element.trim().to_ascii_uppercase().as_str() {
        "H" | "D" => (0.106, 0.0657),
        "C" => (0.339, 0.3598),
        "N" => (0.325, 0.7113),
        "O" => (0.296, 0.8786),
        "F" => (0.312, 0.2552),
        "P" => (0.374, 0.8368),
        "S" => (0.356, 1.0460),
        "CL" => (0.347, 1.1087),
        "NA" => (0.243, 0.0620),
        "MG" => (0.141, 3.7434),
        "K" => (0.303, 0.0137),
        "CA" => (0.243, 1.8800),
        "FE" => (0.227, 0.0556),
        "ZN" => (0.196, 0.0523),
        "BR" => (0.391, 1.2552),
        "I" => (0.418, 1.6736),
        // Carbon-like fallback — keeps the build honest for an exotic
        // element instead of erroring.
        _ => (0.339, 0.3598),
    };
    // The literals above are all valid; `unwrap_or` keeps this
    // total without an `expect` in panel code.
    LjParam::new(sigma, epsilon).unwrap_or(LjParam {
        sigma: 0.339,
        epsilon: 0.3598,
    })
}

/// A built bonded topology — bond / angle / dihedral index lists
/// derived from a bond list. Kept as a plain struct so the assembly is
/// unit-testable without a `System`.
#[derive(Debug, Default, PartialEq)]
struct BondedTopology {
    bonds: Vec<(usize, usize)>,
    angles: Vec<(usize, usize, usize)>,
    dihedrals: Vec<(usize, usize, usize, usize)>,
}

/// Derive angles and proper dihedrals from a bond list.
///
/// - An **angle** `i-j-k` exists for every pair of bonds sharing a
///   central atom `j`.
/// - A proper **dihedral** `i-j-k-l` exists for every bond `j-k` whose
///   two ends each carry at least one further distinct neighbour.
///
/// Indices in each tuple are de-duplicated and ordered canonically so
/// the lists are stable. `O(atoms · degree²)` — fast for the
/// structure sizes a desktop MD panel handles.
fn derive_bonded_topology(n_atoms: usize, bonds: &[(usize, usize)]) -> BondedTopology {
    // Adjacency list.
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
    for &(i, j) in bonds {
        if i < n_atoms && j < n_atoms && i != j {
            if !adj[i].contains(&j) {
                adj[i].push(j);
            }
            if !adj[j].contains(&i) {
                adj[j].push(i);
            }
        }
    }

    // Angles: every (i, k) pair among atom j's neighbours, with i < k.
    let mut angles = Vec::new();
    for (j, nbrs) in adj.iter().enumerate() {
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let (i, k) = (nbrs[a].min(nbrs[b]), nbrs[a].max(nbrs[b]));
                angles.push((i, j, k));
            }
        }
    }

    // Dihedrals: for each bond j-k, pick one i ∈ nbr(j)\{k} and one
    // l ∈ nbr(k)\{j} with i ≠ l. One representative dihedral per
    // central bond keeps the parameter count modest (a real torsion
    // term per central bond is enough for a generic force field).
    let mut dihedrals = Vec::new();
    for &(j, k) in bonds {
        if j >= n_atoms || k >= n_atoms || j == k {
            continue;
        }
        let i = adj[j].iter().copied().find(|&x| x != k);
        let l = adj[k].iter().copied().find(|&x| x != j && Some(x) != i);
        if let (Some(i), Some(l)) = (i, l) {
            dihedrals.push((i, j, k, l));
        }
    }

    BondedTopology {
        bonds: bonds.to_vec(),
        angles,
        dihedrals,
    }
}

/// Build a fully parameterised [`System`] + [`ForceField`] from raw
/// structure text.
///
/// Parses the structure, detects bonds with the covalent-radius rule
/// (scaled by `bond_tolerance`), derives the angle / dihedral
/// topology, writes that connectivity onto the system's
/// [`valenx_md::Topology`], and builds a generic harmonic force field
/// with element-specific Lennard-Jones parameters.
///
/// Returns the system, its force field and the bond / angle /
/// dihedral counts (for the result summary).
fn build_md_system(
    text: &str,
    is_xyz: bool,
    bond_tolerance: f64,
) -> Result<(System, ForceField, [usize; 3]), String> {
    let mut system = if is_xyz {
        read_xyz(text).map_err(|e| format!("parse XYZ: {e}"))?
    } else {
        read_pdb(text).map_err(|e| format!("parse PDB: {e}"))?
    };
    let n = system.len();
    if n < 2 {
        return Err(format!("need at least 2 atoms, got {n}"));
    }

    // --- Bond detection (covalent-radius rule) ---------------------
    // Positions are nm in the System; the molecule_view detector works
    // in ångström — convert via the ViewMolecule bridge.
    let view = ViewMolecule::from_md_system(&system);
    let bonds = detect_bonds_scaled(&view, bond_tolerance);
    if bonds.is_empty() {
        return Err(
            "no bonds detected — atoms are too far apart for the covalent-radius \
             rule (raise the bond tolerance, or check the structure units)"
                .to_string(),
        );
    }
    let bonded = derive_bonded_topology(n, &bonds);

    // --- Write the connectivity onto the topology ------------------
    for &(i, j) in &bonded.bonds {
        system
            .topology
            .add_bond(i, j)
            .map_err(|e| format!("add bond: {e}"))?;
    }
    for &(i, j, k) in &bonded.angles {
        system
            .topology
            .add_angle(i, j, k)
            .map_err(|e| format!("add angle: {e}"))?;
    }
    for &(i, j, k, l) in &bonded.dihedrals {
        system
            .topology
            .add_dihedral(i, j, k, l)
            .map_err(|e| format!("add dihedral: {e}"))?;
    }

    // --- Build the force field -------------------------------------
    let mut ff = ForceField::new(CombiningRule::LorentzBerthelot);
    // Element-specific LJ for every distinct atom type.
    let mut seen: Vec<String> = Vec::new();
    for atom in &system.topology.atoms {
        if !seen.contains(&atom.type_name) {
            seen.push(atom.type_name.clone());
            let element = if atom.element.is_empty() {
                atom.type_name.as_str()
            } else {
                atom.element.as_str()
            };
            ff.set_lj(atom.type_name.clone(), element_lj(element));
        }
    }
    // Generic harmonic bonded parameters. Equilibrium lengths /
    // angles are taken from the *current* geometry so the structure
    // starts near a force-field minimum (a real parameterisation
    // would use type-specific tables; this keeps the v1 honest and
    // stable). Force constants are organic-chemistry-typical values.
    for &(i, j) in &bonded.bonds {
        let r0 = (system.positions[i] - system.positions[j]).norm();
        ff.push_bond(
            BondParam::new(r0.max(0.05), 250_000.0).map_err(|e| format!("bond param: {e}"))?,
        );
    }
    for &(i, j, k) in &bonded.angles {
        let theta0 = bond_angle(
            system.positions[i],
            system.positions[j],
            system.positions[k],
        );
        ff.push_angle(AngleParam::new(theta0, 400.0).map_err(|e| format!("angle param: {e}"))?);
    }
    for _ in &bonded.dihedrals {
        // A mild 3-fold periodic torsion — the generic alkane-like
        // barrier. Phase 0, multiplicity 3, modest height.
        ff.push_dihedral(
            DihedralParam::periodic(2.0, 3, 0.0).map_err(|e| format!("dihedral param: {e}"))?,
        );
    }

    let counts = [
        bonded.bonds.len(),
        bonded.angles.len(),
        bonded.dihedrals.len(),
    ];
    Ok((system, ff, counts))
}

/// Angle (radians) at vertex `b` of the triangle `a-b-c`.
fn bond_angle(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> f64 {
    let u = a - b;
    let v = c - b;
    let nu = u.norm();
    let nv = v.norm();
    if nu < 1e-9 || nv < 1e-9 {
        return std::f64::consts::FRAC_PI_2;
    }
    (u.dot(&v) / (nu * nv)).clamp(-1.0, 1.0).acos()
}

/// Run the bond detector with a tolerance multiplier on the cutoff.
///
/// [`molecule_view::detect_bonds`] uses a fixed slack; here the panel
/// lets the user widen / tighten it by scaling each pair's covalent-
/// radius sum. `tolerance` of `1.0` reproduces the standard rule.
fn detect_bonds_scaled(mol: &ViewMolecule, tolerance: f64) -> Vec<(usize, usize)> {
    let extra = (0.45 * tolerance) as f32;
    let mut bonds = Vec::new();
    for i in 0..mol.atoms.len() {
        let ri = molecule_view::covalent_radius(&mol.atoms[i].element);
        for j in (i + 1)..mol.atoms.len() {
            let rj = molecule_view::covalent_radius(&mol.atoms[j].element);
            let p = mol.atoms[i].pos;
            let q = mol.atoms[j].pos;
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            if d > 0.4 && d <= ri + rj + extra {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Draw Maxwell-Boltzmann velocities into `system` at `temperature`.
///
/// Each velocity component is a Gaussian with standard deviation
/// `√(k_B·T / m)` (the equipartition result). Net centre-of-mass
/// drift is then removed so total momentum starts at zero.
fn seed_maxwell_velocities(system: &mut System, temperature: f64, seed: u64) {
    let mut rng = Rng::new(seed);
    let kb = valenx_md::units::BOLTZMANN;
    let velocities: Vec<Vector3<f64>> = system
        .topology
        .atoms
        .iter()
        .map(|a| {
            let sigma = (kb * temperature / a.mass).max(0.0).sqrt();
            Vector3::new(
                sigma * rng.normal(),
                sigma * rng.normal(),
                sigma * rng.normal(),
            )
        })
        .collect();
    // `set_velocities` only fails on a length mismatch, which the
    // construction above cannot produce.
    let _ = system.set_velocities(velocities);
    system.remove_com_motion();
}

/// Build the chosen integrator for the run.
fn make_integrator(choice: IntegratorChoice, dt_ps: f64) -> Result<Box<dyn Integrator>, String> {
    match choice {
        IntegratorChoice::VelocityVerlet => VelocityVerlet::new(dt_ps)
            .map(|i| Box::new(i) as Box<dyn Integrator>)
            .map_err(|e| e.to_string()),
        IntegratorChoice::Leapfrog => LeapFrog::new(dt_ps)
            .map(|i| Box::new(i) as Box<dyn Integrator>)
            .map_err(|e| e.to_string()),
    }
}

/// Build the chosen thermostat for the run, or `None` for NVE.
fn make_thermostat(
    choice: ThermostatChoice,
    temperature: f64,
    seed: u64,
) -> Result<Option<Box<dyn Thermostat>>, String> {
    match choice {
        ThermostatChoice::None => Ok(None),
        ThermostatChoice::Berendsen => Berendsen::new(temperature, 0.1)
            .map(|t| Some(Box::new(t) as Box<dyn Thermostat>))
            .map_err(|e| e.to_string()),
        ThermostatChoice::Andersen => Andersen::new(temperature, 5.0, seed)
            .map(|t| Some(Box::new(t) as Box<dyn Thermostat>))
            .map_err(|e| e.to_string()),
        ThermostatChoice::VelocityRescale => VelocityRescale::new(temperature, 0.1, seed)
            .map(|t| Some(Box::new(t) as Box<dyn Thermostat>))
            .map_err(|e| e.to_string()),
    }
}

/// Render the Molecular Dynamics panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.md;

    common::section(ui, "Structure");
    ui.horizontal(|ui| {
        ui.radio_value(&mut p.is_xyz, true, "XYZ")
            .on_hover_text("Parse as XYZ — bare element + Cartesian coords, no topology.");
        ui.radio_value(&mut p.is_xyz, false, "PDB")
            .on_hover_text("Parse as PDB — ATOM/HETATM records (CONECT is ignored).");
        if ui
            .small_button("Load file…")
            .on_hover_text("Open a .pdb or .xyz structure file from disk.")
            .clicked()
        {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("Structure", &["pdb", "xyz"])
                .pick_file()
            {
                // Round-21 H1: see biostruct loader.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => {
                        p.is_xyz = path
                            .extension()
                            .map(|e| e.eq_ignore_ascii_case("xyz"))
                            .unwrap_or(p.is_xyz);
                        p.structure_text = t;
                    }
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
    });
    ui.add(
        egui::TextEdit::multiline(&mut p.structure_text)
            .id_source("md_structure_input")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(5),
    );
    ui.horizontal(|ui| {
        let lbl = ui.label("Bond-detection tolerance:").on_hover_text(
            "Multiplier on the covalent-radius cutoff used to infer bonds from coordinates.",
        );
        ui.add(
            egui::DragValue::new(&mut p.bond_tolerance)
                .speed(0.05)
                .range(0.6..=1.6),
        )
        .labelled_by(lbl.id)
        .on_hover_text(
            "1.0 = textbook covalent-radius rule. Bump to 1.2 if your input has noisy \
             coordinates; drop to 0.9 to be stricter about long bonds.",
        );
        ui.label("× covalent radii");
    });
    ui.separator();

    common::section(ui, "Force field");
    ui.label(
        egui::RichText::new(
            "bonds + angles + dihedrals from detected topology · \
             element-specific Lennard-Jones",
        )
        .weak()
        .small(),
    );
    ui.separator();

    common::section(ui, "Run configuration");
    egui::Grid::new("md_params")
        .num_columns(2)
        .spacing([8.0, 4.0])
        .show(ui, |ui| {
            ui.label("Integrator");
            egui::ComboBox::from_id_source("md_integrator")
                .selected_text(p.integrator.label())
                .show_ui(ui, |ui| {
                    for c in [IntegratorChoice::VelocityVerlet, IntegratorChoice::Leapfrog] {
                        ui.selectable_value(&mut p.integrator, c, c.label());
                    }
                });
            ui.end_row();

            ui.label("Thermostat");
            egui::ComboBox::from_id_source("md_thermostat")
                .selected_text(p.thermostat.label())
                .show_ui(ui, |ui| {
                    for c in [
                        ThermostatChoice::None,
                        ThermostatChoice::Berendsen,
                        ThermostatChoice::Andersen,
                        ThermostatChoice::VelocityRescale,
                    ] {
                        ui.selectable_value(&mut p.thermostat, c, c.label());
                    }
                });
            ui.end_row();

            let lbl = ui
                .label("Temperature (K)")
                .on_hover_text("Thermostat setpoint. SI unit: K. 300 K ≈ room temperature.");
            ui.add(
                egui::DragValue::new(&mut p.temperature)
                    .speed(5.0)
                    .range(1.0..=1000.0),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Thermostat target temperature in Kelvin. Initial velocities are \
                 drawn from Maxwell-Boltzmann at this T.",
            );
            ui.end_row();

            let lbl = ui
                .label("Time step (fs)")
                .on_hover_text("Integration step. SI unit: fs (10⁻¹⁵ s).");
            ui.add(
                egui::DragValue::new(&mut p.timestep_fs)
                    .speed(0.1)
                    .range(0.1..=4.0),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Integration time step in femtoseconds. 1.0 fs is the textbook all-atom \
                 default; ≤ 2.0 fs without hydrogen constraints; ≤ 4.0 fs with SHAKE.",
            );
            ui.end_row();

            let lbl = ui
                .label("Steps")
                .on_hover_text("Total number of integration steps.");
            ui.add(egui::DragValue::new(&mut p.steps).range(1..=20_000))
                .labelled_by(lbl.id)
                .on_hover_text("Number of MD steps. Total simulated time = steps × time-step.");
            ui.end_row();

            let lbl = ui
                .label("Report every")
                .on_hover_text("Energy / temperature sampling interval (in steps).");
            ui.add(egui::DragValue::new(&mut p.report_interval).range(1..=1000))
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Number of steps between energy / temperature samples in the report.",
                );
            ui.end_row();

            let lbl = ui.label("RNG seed").on_hover_text(
                "Pseudorandom seed for initial-velocity sampling + stochastic thermostats. \
                     Same seed = exactly reproducible run.",
            );
            ui.add(egui::DragValue::new(&mut p.seed))
                .labelled_by(lbl.id)
                .on_hover_text("Reproducibility seed.");
            ui.end_row();
        });

    if common::run_button(ui, "Run simulation") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_simulation(p);
    }
    ui.horizontal(|ui| {
        let (u, r) = common::undo_redo_inline(ui, p.can_undo(), p.can_redo());
        if u {
            p.undo_edit();
        }
        if r {
            p.redo_edit();
        }
        ui.label(
            egui::RichText::new("Ctrl+Z / Ctrl+Y reverses last Run")
                .weak()
                .small(),
        );
    });

    common::error_line(ui, &p.error);

    // --- Push the built structure into the 3-D viewport -----------
    if !p.structure_text.trim().is_empty() {
        ui.horizontal(|ui| {
            if ui.button("Show in 3D viewport").clicked() {
                show_in_viewport(app, false);
            }
            if ui.button("Show (spacefill)").clicked() {
                show_in_viewport(app, true);
            }
        });
    }
    // `app.genetics.md` is re-borrowed after the viewport calls above.
    let p = &mut app.genetics.md;

    if !p.trace.is_empty() {
        ui.separator();
        common::section(ui, "Energy trace (kJ/mol)");
        let pot: PlotPoints = p.trace.iter().map(|(s, v, _, _)| [*s, *v]).collect();
        let kin: PlotPoints = p.trace.iter().map(|(s, _, v, _)| [*s, *v]).collect();
        let tot: PlotPoints = p.trace.iter().map(|(s, _, _, v)| [*s, *v]).collect();
        managed_plot_mem_cfg(
            ui,
            "md_energy_plot",
            150.0,
            |plot| plot.legend(Legend::default()),
            |plot_ui| {
                plot_ui.line(Line::new(pot).name("potential"));
                plot_ui.line(Line::new(kin).name("kinetic"));
                plot_ui.line(Line::new(tot).name("total"));
            },
        );
    }

    if !p.temp_trace.is_empty() {
        common::section(ui, "Temperature (K)");
        let temp: PlotPoints = p.temp_trace.iter().map(|(s, t)| [*s, *t]).collect();
        managed_plot_mem_cfg(
            ui,
            "md_temp_plot",
            120.0,
            |plot| plot.legend(Legend::default()),
            |plot_ui| {
                plot_ui.line(Line::new(temp).name("temperature"));
            },
        );
    }

    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "md_result", &p.result, 13);
    }
}

/// Execute the configured MD run, filling the panel's traces + result.
fn run_simulation(p: &mut MdPanel) {
    p.error = None;
    p.trace.clear();
    p.temp_trace.clear();

    let (mut system, ff, counts) =
        match build_md_system(&p.structure_text, p.is_xyz, p.bond_tolerance) {
            Ok(v) => v,
            Err(e) => {
                p.error = Some(e);
                return;
            }
        };
    let n_atoms = system.len();
    let start_positions = system.positions.clone();

    // Maxwell-Boltzmann initial velocities at the target temperature.
    seed_maxwell_velocities(&mut system, p.temperature, p.seed);

    let dt_ps = p.timestep_fs * valenx_md::units::FS_TO_PS;
    let integrator = match make_integrator(p.integrator, dt_ps) {
        Ok(i) => i,
        Err(e) => {
            p.error = Some(format!("integrator: {e}"));
            return;
        }
    };
    let thermostat = match make_thermostat(p.thermostat, p.temperature, p.seed) {
        Ok(t) => t,
        Err(e) => {
            p.error = Some(format!("thermostat: {e}"));
            return;
        }
    };

    let mut sim = match Simulation::new(system, ff) {
        Ok(s) => s,
        Err(e) => {
            p.error = Some(format!("simulation setup: {e}"));
            return;
        }
    };
    sim = sim.with_integrator(integrator);
    if let Some(t) = thermostat {
        sim = sim.with_thermostat(t);
    }
    if let Err(e) = sim.set_report_interval(p.report_interval) {
        p.error = Some(e.to_string());
        return;
    }

    let report = match sim.run(p.steps) {
        Ok(r) => r,
        Err(e) => {
            p.error = Some(format!("run: {e}"));
            return;
        }
    };

    for r in &sim.log.reports {
        p.trace.push((
            r.step as f64,
            r.potential_energy,
            r.kinetic_energy,
            r.total_energy,
        ));
        p.temp_trace.push((r.step as f64, r.temperature));
    }

    // All-atom RMSD from the starting structure (Kabsch-superposed).
    let final_rmsd = kabsch_rmsd(&sim.system.positions, &start_positions)
        .map(|v| format!("{v:.4} nm"))
        .unwrap_or_else(|e| format!("(n/a: {e})"));

    p.result = format!(
        "atoms          : {n_atoms}\n\
         bonds          : {}\nangles         : {}\ndihedrals      : {}\n\
         integrator     : {}\nthermostat     : {}\n\
         steps          : {}\ntime step      : {:.2} fs\n\
         final time     : {:.4} ps\n\
         potential E    : {:.4} kJ/mol\n\
         kinetic E      : {:.4} kJ/mol\n\
         total E        : {:.4} kJ/mol\n\
         temperature    : {:.2} K  (target {:.0} K)\n\
         energy drift σ : {:.4} kJ/mol\n\
         final RMSD     : {final_rmsd}",
        counts[0],
        counts[1],
        counts[2],
        p.integrator.label(),
        p.thermostat.label(),
        report.steps,
        p.timestep_fs,
        report.final_time,
        report.final_potential_energy,
        report.final_kinetic_energy,
        report.final_total_energy,
        report.final_temperature,
        p.temperature,
        sim.log.total_energy_std(),
    );
}

/// Build the loaded structure's ball-and-stick (or spacefill) mesh and
/// push it into the app's 3-D viewport.
fn show_in_viewport(app: &mut ValenxApp, spacefill: bool) {
    let p = &app.genetics.md;
    let system = if p.is_xyz {
        read_xyz(&p.structure_text)
    } else {
        read_pdb(&p.structure_text)
    };
    let system = match system {
        Ok(s) => s,
        Err(e) => {
            app.genetics.md.error = Some(format!("parse: {e}"));
            return;
        }
    };
    let mut view = ViewMolecule::from_md_system(&system);
    // Honour the panel's bond tolerance so the sticks match the run.
    view.bonds = detect_bonds_scaled(&view, app.genetics.md.bond_tolerance);
    let mesh = if spacefill {
        molecule_view::spacefill(&view)
    } else {
        molecule_view::ball_and_stick(&view, 0.25, 0.16)
    };
    let label = if spacefill {
        "md-structure.spacefill"
    } else {
        "md-structure.ball-stick"
    };
    match molecule_view::show_molecule(app, mesh, label) {
        Ok(_) => app.genetics.md.error = None,
        Err(e) => app.genetics.md.error = Some(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_ethane_builds_a_bonded_system() {
        let p = MdPanel::default();
        let (system, ff, counts) =
            build_md_system(&p.structure_text, p.is_xyz, p.bond_tolerance).expect("build");
        assert_eq!(system.len(), 8);
        // Ethane: 1 C-C + 6 C-H = 7 bonds.
        assert_eq!(counts[0], 7);
        // Angles and dihedrals must both be non-empty for a real
        // bonded MD setup.
        assert!(counts[1] > 0, "no angles derived");
        assert!(counts[2] > 0, "no dihedrals derived");
        // The force field must cover the topology.
        assert!(ff.validate_against(&system.topology).is_ok());
    }

    #[test]
    fn default_run_produces_traces_and_rmsd() {
        let mut p = MdPanel {
            steps: 30,
            report_interval: 5,
            ..MdPanel::default()
        };
        run_simulation(&mut p);
        assert!(p.error.is_none(), "run errored: {:?}", p.error);
        assert!(!p.trace.is_empty(), "energy trace empty");
        assert!(!p.temp_trace.is_empty(), "temperature trace empty");
        assert!(p.result.contains("final RMSD"));
    }

    #[test]
    fn derive_bonded_topology_for_a_methane_star() {
        // C at 0 bonded to four H (1..=4): 4 bonds, C(4,2)=6 angles,
        // no proper dihedral (every dihedral needs a second central
        // atom with its own neighbour).
        let bonds = vec![(0, 1), (0, 2), (0, 3), (0, 4)];
        let bt = derive_bonded_topology(5, &bonds);
        assert_eq!(bt.bonds.len(), 4);
        assert_eq!(bt.angles.len(), 6);
        assert!(bt.dihedrals.is_empty());
    }

    #[test]
    fn derive_bonded_topology_for_a_four_atom_chain() {
        // 0-1-2-3 linear chain: 3 bonds, 2 angles, 1 dihedral.
        let bonds = vec![(0, 1), (1, 2), (2, 3)];
        let bt = derive_bonded_topology(4, &bonds);
        assert_eq!(bt.bonds.len(), 3);
        assert_eq!(bt.angles.len(), 2);
        assert_eq!(bt.dihedrals.len(), 1);
        assert_eq!(bt.dihedrals[0], (0, 1, 2, 3));
    }

    #[test]
    fn element_lj_is_element_specific() {
        // Distinct elements get distinct parameters — not one global
        // σ/ε.
        let c = element_lj("C");
        let h = element_lj("H");
        let o = element_lj("O");
        assert_ne!(c, h);
        assert_ne!(c, o);
        // Hydrogen has the smallest σ of the three.
        assert!(h.sigma < c.sigma && h.sigma < o.sigma);
    }

    #[test]
    fn build_fails_when_atoms_are_unbonded() {
        // Two carbons 10 Å apart — no bond, no MD setup.
        let xyz = "2\nfar\nC 0.0 0.0 0.0\nC 10.0 0.0 0.0\n";
        let err = build_md_system(xyz, true, 1.0).unwrap_err();
        assert!(err.contains("no bonds"), "got: {err}");
    }

    #[test]
    fn seed_maxwell_velocities_warms_a_cold_system() {
        let (mut system, _ff, _c) =
            build_md_system(&MdPanel::default().structure_text, true, 1.0).unwrap();
        assert!(system.temperature(0) < 1e-9, "should start cold");
        seed_maxwell_velocities(&mut system, 300.0, 7);
        // Temperature is now non-trivial and momentum was zeroed.
        assert!(system.temperature(0) > 1.0);
        assert!(system.linear_momentum().norm() < 1e-6);
    }

    #[test]
    fn bond_angle_of_a_right_angle_is_ninety_degrees() {
        let a = Vector3::new(1.0, 0.0, 0.0);
        let b = Vector3::new(0.0, 0.0, 0.0);
        let c = Vector3::new(0.0, 1.0, 0.0);
        assert!((bond_angle(a, b, c) - std::f64::consts::FRAC_PI_2).abs() < 1e-9);
    }
}

/// Headless egui UI-logic tests for the Molecular Dynamics panel.
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

    fn app_with_panel() -> ValenxApp {
        let mut app = ValenxApp::default();
        app.genetics.active = GeneticsPanel::MolecularDynamics;
        app
    }

    #[test]
    fn draws_fresh_state_without_panic() {
        let mut app = app_with_panel();
        draw_headless(&mut app);
    }

    #[test]
    fn draws_post_run_state_with_traces_without_panic() {
        // Post-run: energy + temperature traces drive egui_plot line
        // charts; a result string drives the mono output.
        let mut app = app_with_panel();
        app.genetics.md.trace = vec![(0.0, -1.0, 0.5, -0.5), (1.0, -1.1, 0.6, -0.5)];
        app.genetics.md.temp_trace = vec![(0.0, 280.0), (1.0, 305.0)];
        app.genetics.md.result = "atoms : 8\nfinal RMSD : 0.0123 nm\n".to_string();
        draw_headless(&mut app);
        // Error state.
        let mut app = app_with_panel();
        app.genetics.md.error = Some("no bonds detected".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_simulation_produces_traces_and_a_result() {
        // A short MD run calls the real valenx-md API and fills the
        // energy / temperature traces + result.
        let mut p = MdPanel {
            steps: 20,
            report_interval: 5,
            ..MdPanel::default()
        };
        run_simulation(&mut p);
        assert!(p.error.is_none(), "MD run errored: {:?}", p.error);
        assert!(!p.trace.is_empty(), "energy trace empty");
        assert!(!p.temp_trace.is_empty(), "temperature trace empty");
        assert!(p.result.contains("final RMSD"));
        assert!(p.result.contains("integrator"));
    }

    #[test]
    fn run_simulation_surfaces_error_on_unbonded_structure() {
        // Two carbons 10 Å apart cannot form a bonded MD system — the
        // panel must surface an error rather than panicking.
        let mut p = MdPanel {
            structure_text: "2\nfar\nC 0.0 0.0 0.0\nC 10.0 0.0 0.0\n".to_string(),
            is_xyz: true,
            steps: 10,
            ..MdPanel::default()
        };
        run_simulation(&mut p);
        assert!(
            p.error.is_some(),
            "MD run should error on an unbonded structure"
        );
        assert!(p.trace.is_empty());
    }

    #[test]
    fn run_simulation_surfaces_error_on_malformed_structure() {
        // Junk text is not a valid XYZ structure.
        let mut p = MdPanel {
            structure_text: "this is not a structure file".to_string(),
            is_xyz: true,
            steps: 10,
            ..MdPanel::default()
        };
        run_simulation(&mut p);
        assert!(p.error.is_some(), "MD run should error on malformed input");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        use egui::accesskit::Role;
        let mut app = app_with_panel();
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                super::draw(&mut app, ui);
            });
        });
        let nodes = out
            .platform_output
            .accesskit_update
            .expect("accesskit tree produced")
            .nodes;
        let spin_buttons: Vec<_> = nodes
            .iter()
            .filter(|(_, n)| n.role() == Role::SpinButton)
            .collect();
        assert!(
            !spin_buttons.is_empty(),
            "MD panel should expose at least one SpinButton"
        );
        assert!(
            spin_buttons
                .iter()
                .all(|(_, n)| !n.labelled_by().is_empty()),
            "every MD DragValue must be labelled_by its caption (AI-drivable name)"
        );
    }
}
