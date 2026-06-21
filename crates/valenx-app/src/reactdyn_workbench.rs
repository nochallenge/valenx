//! The right-side **Reaction Dynamics** workbench — native 3-D ab-initio
//! molecular dynamics over `valenx-reactdyn` (Born-Oppenheimer AIMD:
//! velocity-Verlet + numerical qchem forces).
//!
//! Pick a preset reaction or paste an XYZ, choose the method / basis /
//! timestep / steps / thermostat, and run. The solve is small-system +
//! numerical-gradient, so it runs on a **background thread** (the UI
//! stays responsive) and is **cost-guarded** (atoms × steps capped) so it
//! can never lock up the machine. The finished trajectory plays back in
//! the **3-D viewport** — atoms move and **bonds form/break** (recomputed
//! each frame) — alongside an energy plot whose flat total-energy line is
//! the honest correctness check.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use valenx_qchem::element::Element;
use valenx_qchem::geometry::{MolecularGeometry, BOHR_PER_ANGSTROM};
use valenx_reactdyn::{
    morse_param, AimdEngine, Controls, Embedding, Method, MmAtom, QmMmEngine, QmMmSystem,
    ReactionEngine, ReactiveEngine, ReactiveSystem, System, Thermostat, Trajectory,
};

use crate::genetics::molecule_view::{self, ViewAtom, ViewMolecule};
use crate::ValenxApp;

/// A built-in small-molecule starting point.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Preset {
    /// H₂ started a little stretched — vibrates around the bond.
    #[default]
    H2,
    /// Water (H₂O) — a 3-atom flex demo.
    Water,
    /// Hydrogen fluoride (HF) — a polar diatomic.
    Hf,
    /// User-supplied XYZ.
    Custom,
}

impl Preset {
    const ALL: [Preset; 4] = [Preset::H2, Preset::Water, Preset::Hf, Preset::Custom];

    fn label(self) -> &'static str {
        match self {
            Preset::H2 => "H₂ (stretched)",
            Preset::Water => "Water (H₂O)",
            Preset::Hf => "Hydrogen fluoride (HF)",
            Preset::Custom => "Custom (XYZ)",
        }
    }

    /// The starting geometry as XYZ text (ångström). `None` for Custom.
    fn xyz(self) -> Option<&'static str> {
        match self {
            Preset::H2 => Some("2\nH2 stretched\nH 0.0 0.0 0.0\nH 0.0 0.0 0.95\n"),
            Preset::Water => Some(
                "3\nwater  charge 0  mult 1\n\
                 O  0.000000  0.000000  0.117300\n\
                 H  0.000000  0.757200 -0.469200\n\
                 H  0.000000 -0.757200 -0.469200\n",
            ),
            Preset::Hf => Some("2\nHF\nH 0.0 0.0 0.0\nF 0.0 0.0 0.95\n"),
            Preset::Custom => None,
        }
    }
}

/// Which thermostat the form has selected.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum ThermostatKind {
    /// Microcanonical (energy-conserving).
    #[default]
    Nve,
    /// Berendsen weak coupling to a target temperature.
    Berendsen,
}

/// Which physics backend the run uses.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
enum Backend {
    /// Ab-initio MD — the molecule alone, in vacuum.
    #[default]
    Aimd,
    /// QM/MM — the molecule (QM) in an explicit classical solvent shell.
    QmMm,
    /// Reactive classical force field — a many-atom cluster (materials).
    Reactive,
}

/// A live background AIMD run.
struct RunHandle {
    /// Filled by the worker thread when the run finishes (Ok or Err).
    result: Arc<Mutex<Option<Result<Trajectory, String>>>>,
    /// Steps completed so far (for the progress bar).
    progress: Arc<AtomicUsize>,
    /// Total steps requested.
    total: usize,
    /// Keeps the worker joined when dropped.
    _handle: thread::JoinHandle<()>,
}

/// Form + result state for the Reaction Dynamics workbench.
pub struct ReactdynWorkbenchState {
    backend: Backend,
    preset: Preset,
    xyz: String,
    /// Number of explicit MM solvent atoms (QM/MM backend).
    n_solvent: usize,
    /// Number of atoms in the reactive cluster (Reactive backend).
    n_cluster: usize,
    /// QM/MM embedding scheme.
    embedding: Embedding,
    method: Method,
    basis: String,
    dt_fs: f64,
    n_steps: usize,
    thermostat: ThermostatKind,
    target_kelvin: f64,
    tau_fs: f64,
    run: Option<RunHandle>,
    last: Option<Trajectory>,
    status: String,
    error: Option<String>,
    // Playback.
    frame_idx: usize,
    playing: bool,
    /// The frame index currently pushed to the viewport (avoids
    /// rebuilding/re-pushing the mesh on every repaint).
    last_pushed: Option<usize>,
}

impl Default for ReactdynWorkbenchState {
    fn default() -> Self {
        Self {
            backend: Backend::Aimd,
            preset: Preset::H2,
            xyz: String::new(),
            n_solvent: 10,
            n_cluster: 12,
            embedding: Embedding::Mechanical,
            method: Method::Rhf,
            basis: "STO-3G".to_string(),
            dt_fs: 0.5,
            n_steps: 40,
            thermostat: ThermostatKind::Nve,
            target_kelvin: 300.0,
            tau_fs: 50.0,
            run: None,
            last: None,
            status: String::new(),
            error: None,
            frame_idx: 0,
            playing: false,
            last_pushed: None,
        }
    }
}

impl ReactdynWorkbenchState {
    /// `true` while a background run is in flight.
    pub fn is_running(&self) -> bool {
        self.run.is_some()
    }

    fn build_thermostat(&self) -> Thermostat {
        match self.thermostat {
            ThermostatKind::Nve => Thermostat::Nve,
            ThermostatKind::Berendsen => Thermostat::Berendsen {
                target_kelvin: self.target_kelvin,
                tau_fs: self.tau_fs,
            },
        }
    }
}

/// Parse the selected setup (preset or custom XYZ) into a geometry.
fn parse_geometry(s: &ReactdynWorkbenchState) -> Result<MolecularGeometry, String> {
    let xyz_owned;
    let xyz: &str = match s.preset.xyz() {
        Some(text) => text,
        None => {
            xyz_owned = s.xyz.clone();
            xyz_owned.as_str()
        }
    };
    let geom = MolecularGeometry::from_xyz_str(xyz).map_err(|e| format!("geometry: {e}"))?;
    if geom.atoms.is_empty() {
        return Err("geometry has no atoms".into());
    }
    Ok(geom)
}

/// The shared, cost-guarded run controls.
fn build_controls(s: &ReactdynWorkbenchState) -> Controls {
    Controls {
        method: s.method,
        basis: s.basis.clone(),
        dt_fs: s.dt_fs,
        n_steps: s.n_steps,
        fd_delta_bohr: 0.01,
        thermostat: s.build_thermostat(),
        // Cap QM atoms × steps so a run is always bounded and safe.
        max_cost_guard: 6000,
    }
}

/// Build the AIMD (vacuum) inputs. Extracted for headless tests.
fn build_aimd_inputs(s: &ReactdynWorkbenchState) -> Result<(System, Controls), String> {
    let geom = parse_geometry(s)?;
    let system = System {
        elements: geom.atoms.iter().map(|a| a.element).collect(),
        pos_bohr: geom.atoms.iter().map(|a| a.position).collect(),
        charge: geom.charge,
        multiplicity: geom.multiplicity,
    };
    Ok((system, build_controls(s)))
}

/// Build the QM/MM inputs: the molecule as the QM region + an explicit
/// Lennard-Jones (neon) solvent shell placed evenly around it.
fn build_qmmm_inputs(s: &ReactdynWorkbenchState) -> Result<(QmMmSystem, Controls), String> {
    let geom = parse_geometry(s)?;
    let qm_elements: Vec<Element> = geom.atoms.iter().map(|a| a.element).collect();
    let qm_pos_bohr: Vec<[f64; 3]> = geom.atoms.iter().map(|a| a.position).collect();
    let qm_classical: Vec<(f64, f64, f64)> = qm_elements
        .iter()
        .map(|e| {
            let (sigma, epsilon) = lj_params(e.symbol());
            (0.0, sigma, epsilon) // mechanical LJ coupling; QM charge 0 for v1
        })
        .collect();

    // Centroid of the QM region (bohr).
    let nq = qm_pos_bohr.len() as f64;
    let mut centroid = [0.0_f64; 3];
    for p in &qm_pos_bohr {
        centroid[0] += p[0];
        centroid[1] += p[1];
        centroid[2] += p[2];
    }
    centroid[0] /= nq;
    centroid[1] /= nq;
    centroid[2] /= nq;

    // A neon LJ solvent shell ~5 Å out, evenly spread (no RNG).
    let ne = Element::from_symbol("Ne").map_err(|e| format!("element: {e}"))?;
    let (ne_sigma, ne_eps) = lj_params("Ne");
    let radius_bohr = 5.0 * BOHR_PER_ANGSTROM;
    let mm: Vec<MmAtom> = fibonacci_sphere(s.n_solvent, radius_bohr, centroid)
        .into_iter()
        .enumerate()
        .map(|(i, pos_bohr)| MmAtom {
            element: ne,
            pos_bohr,
            // A model polar solvent: alternating +/- partial charges
            // (net ~0). Only affects the electrostatic-embedding path.
            charge: if i % 2 == 0 { 0.3 } else { -0.3 },
            sigma_bohr: ne_sigma,
            epsilon_hartree: ne_eps,
            mass_amu: ne.atomic_mass(),
        })
        .collect();

    let qmmm = QmMmSystem {
        qm_elements,
        qm_pos_bohr,
        qm_charge: geom.charge,
        qm_mult: geom.multiplicity,
        qm_classical,
        mm,
        embedding: s.embedding,
    };
    Ok((qmmm, build_controls(s)))
}

/// Build a reactive-cluster run: a carbon ring driven by the Morse
/// reactive force field. The cost guard is generous — the potential is
/// classical (cheap), so many atoms × steps are fine.
fn build_reactive_inputs(s: &ReactdynWorkbenchState) -> Result<(ReactiveSystem, Controls), String> {
    let n = s.n_cluster.max(2);
    let c = Element::from_symbol("C").map_err(|e| format!("element: {e}"))?;
    let re = morse_param("C").r_e; // bohr
                                   // A planar ring with adjacent atoms one equilibrium bond apart.
    let radius = re / (2.0 * (std::f64::consts::PI / n as f64).sin());
    let pos_bohr: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let a = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            [radius * a.cos(), radius * a.sin(), 0.0]
        })
        .collect();
    let controls = Controls {
        max_cost_guard: 200_000, // classical: cheap, allow many atoms × steps
        ..build_controls(s)
    };
    Ok((
        ReactiveSystem {
            elements: vec![c; n],
            pos_bohr,
        },
        controls,
    ))
}

/// Approximate Lennard-Jones parameters `(σ bohr, ε hartree)` by element
/// — UFF-class values for the mechanical coupling. A coarse default
/// covers anything not listed.
fn lj_params(symbol: &str) -> (f64, f64) {
    // (σ in Å, ε in kcal/mol).
    let (sigma_a, eps_kcal) = match symbol {
        "H" => (2.50, 0.030),
        "C" => (3.40, 0.086),
        "N" => (3.25, 0.170),
        "O" => (3.12, 0.210),
        "F" => (3.00, 0.061),
        "Ne" => (2.78, 0.069),
        _ => (3.00, 0.050),
    };
    const KCAL_MOL_TO_HARTREE: f64 = 0.001_593_601;
    (sigma_a * BOHR_PER_ANGSTROM, eps_kcal * KCAL_MOL_TO_HARTREE)
}

/// Place `n` points evenly on a sphere of `radius` about `center` via the
/// Fibonacci-spiral construction (fully deterministic — no RNG).
fn fibonacci_sphere(n: usize, radius: f64, center: [f64; 3]) -> Vec<[f64; 3]> {
    if n == 0 {
        return Vec::new();
    }
    let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..n)
        .map(|i| {
            let y = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = (1.0 - y * y).max(0.0).sqrt();
            let theta = golden * i as f64;
            [
                center[0] + radius * r * theta.cos(),
                center[1] + radius * y,
                center[2] + radius * r * theta.sin(),
            ]
        })
        .collect()
}

/// Spawn the background run for the selected backend (AIMD or QM/MM).
fn start_run(s: &mut ReactdynWorkbenchState) {
    s.error = None;

    enum Job {
        Aimd(System, Controls),
        QmMm(QmMmSystem, Controls),
        Reactive(ReactiveSystem, Controls),
    }
    let job = match s.backend {
        Backend::Aimd => match build_aimd_inputs(s) {
            Ok((sys, c)) => Job::Aimd(sys, c),
            Err(e) => {
                s.error = Some(e);
                return;
            }
        },
        Backend::QmMm => match build_qmmm_inputs(s) {
            Ok((sys, c)) => Job::QmMm(sys, c),
            Err(e) => {
                s.error = Some(e);
                return;
            }
        },
        Backend::Reactive => match build_reactive_inputs(s) {
            Ok((sys, c)) => Job::Reactive(sys, c),
            Err(e) => {
                s.error = Some(e);
                return;
            }
        },
    };
    let total = match &job {
        Job::Aimd(_, c) | Job::QmMm(_, c) | Job::Reactive(_, c) => c.n_steps,
    };

    let result = Arc::new(Mutex::new(None));
    let progress = Arc::new(AtomicUsize::new(0));
    let (result_w, progress_w) = (Arc::clone(&result), Arc::clone(&progress));
    let handle = thread::spawn(move || {
        let mut on_step = |step: usize| progress_w.store(step + 1, Ordering::Relaxed);
        let outcome = match job {
            Job::Aimd(sys, c) => AimdEngine
                .run(&sys, &c, &mut on_step)
                .map_err(|e| e.to_string()),
            Job::QmMm(sys, c) => QmMmEngine
                .run(&sys, &c, &mut on_step)
                .map_err(|e| e.to_string()),
            Job::Reactive(sys, c) => ReactiveEngine
                .run(&sys, &c, &mut on_step)
                .map_err(|e| e.to_string()),
        };
        *result_w.lock().unwrap() = Some(outcome);
    });
    s.status = format!("running {total} steps…");
    s.run = Some(RunHandle {
        result,
        progress,
        total,
        _handle: handle,
    });
}

/// Move a finished background run's result into `last` / `error`.
fn poll_run(s: &mut ReactdynWorkbenchState) {
    let done = if let Some(run) = &s.run {
        run.result.lock().unwrap().is_some()
    } else {
        false
    };
    if done {
        let run = s.run.take().unwrap();
        let outcome = run.result.lock().unwrap().take();
        match outcome {
            Some(Ok(traj)) => {
                s.status = format!("done — {} frames", traj.frames.len());
                s.last = Some(traj);
                // Reset playback to the start of the new trajectory.
                s.frame_idx = 0;
                s.playing = false;
                s.last_pushed = None;
            }
            Some(Err(e)) => {
                s.status = "failed".into();
                s.error = Some(e);
            }
            None => {}
        }
    }
}

/// Draw the Reaction Dynamics workbench right-side panel.
pub fn draw_reactdyn_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_reactdyn_workbench {
        return;
    }
    poll_run(&mut app.reactdyn);

    // A mesh to push into the 3-D viewport, set inside the panel closure
    // and applied *after* it releases its &mut borrow of `app`.
    let mut to_push: Option<(ViewMolecule, String)> = None;

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_reactdyn_workbench",
        "Reaction Dynamics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native ab-initio MD · valenx-reactdyn")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.reactdyn;
            let running = s.run.is_some();

            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("System").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.backend, Backend::Aimd, "AIMD (vacuum)")
                            .on_hover_text("The molecule alone, in vacuum.");
                        ui.radio_value(&mut s.backend, Backend::QmMm, "QM/MM (solvent)")
                            .on_hover_text("The molecule (QM) in an explicit classical solvent shell.");
                        ui.radio_value(&mut s.backend, Backend::Reactive, "Reactive (materials)")
                            .on_hover_text("Many-atom classical reactive force field (Morse) — bonds form/break, fast.");
                    });
                    if s.backend == Backend::QmMm {
                        ui.horizontal(|ui| {
                            ui.label("solvent atoms");
                            ui.add(egui::DragValue::new(&mut s.n_solvent).speed(0.5));
                        });
                        ui.horizontal(|ui| {
                            ui.radio_value(&mut s.embedding, Embedding::Mechanical, "mechanical")
                                .on_hover_text("Classical LJ + Coulomb coupling; the solvent does not polarize the QM density.");
                            if ui
                                .radio_value(&mut s.embedding, Embedding::Electrostatic, "electrostatic")
                                .on_hover_text("MM charges enter the SCF and polarize the QM density (RHF only). More accurate, slower.")
                                .clicked()
                            {
                                s.method = Method::Rhf; // electrostatic embedding is RHF-only in v1
                            }
                        });
                        ui.label(
                            egui::RichText::new(
                                "A model polar solvent shell (LJ + alternating ± charges) around \
                                 the molecule. Mechanical = LJ coupling; electrostatic = the MM \
                                 charges polarize the quantum density (RHF, slower).",
                            )
                            .weak()
                            .small(),
                        );
                    }
                    if s.backend == Backend::Reactive {
                        ui.horizontal(|ui| {
                            ui.label("cluster atoms (carbon)");
                            ui.add(egui::DragValue::new(&mut s.n_cluster).speed(0.5));
                        });
                        ui.label(
                            egui::RichText::new(
                                "A carbon ring driven by the reactive Morse force field \
                                 (classical, fast). Bonds form/break by distance — run with the \
                                 Berendsen thermostat at high T to drive rearrangement. The \
                                 method / basis below are ignored (the force field is classical).",
                            )
                            .weak()
                            .small(),
                        );
                    } else {
                        egui::ComboBox::from_id_source("reactdyn_preset")
                            .selected_text(s.preset.label())
                            .show_ui(ui, |ui| {
                                for p in Preset::ALL {
                                    ui.selectable_value(&mut s.preset, p, p.label());
                                }
                            });
                        if s.preset == Preset::Custom {
                            ui.add(
                                egui::TextEdit::multiline(&mut s.xyz)
                                    .id_source("reactdyn_xyz")
                                    .font(egui::TextStyle::Monospace)
                                    .desired_rows(6)
                                    .desired_width(f32::INFINITY)
                                    .hint_text("paste an XYZ geometry (ångström)"),
                            );
                        }
                    }

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Method").strong());
                    ui.horizontal_wrapped(|ui| {
                        ui.radio_value(&mut s.method, Method::Rhf, "RHF");
                        ui.radio_value(&mut s.method, Method::Uhf, "UHF");
                        ui.radio_value(&mut s.method, Method::Dft, "DFT (B3LYP)");
                    });
                    ui.horizontal(|ui| {
                        ui.label("basis");
                        egui::ComboBox::from_id_source("reactdyn_basis")
                            .selected_text(&s.basis)
                            .show_ui(ui, |ui| {
                                for b in ["STO-3G", "3-21G", "6-31G"] {
                                    ui.selectable_value(&mut s.basis, b.to_string(), b);
                                }
                            });
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Dynamics").strong());
                    ui.horizontal(|ui| {
                        ui.label("dt (fs)");
                        ui.add(egui::DragValue::new(&mut s.dt_fs).speed(0.05));
                        ui.label("steps");
                        ui.add(egui::DragValue::new(&mut s.n_steps).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.radio_value(&mut s.thermostat, ThermostatKind::Nve, "NVE")
                            .on_hover_text("Energy-conserving — total energy stays flat (the correctness check).");
                        ui.radio_value(&mut s.thermostat, ThermostatKind::Berendsen, "Berendsen")
                            .on_hover_text("Weak coupling to a target temperature.");
                    });
                    if s.thermostat == ThermostatKind::Berendsen {
                        ui.horizontal(|ui| {
                            ui.label("T (K)");
                            ui.add(egui::DragValue::new(&mut s.target_kelvin).speed(5.0));
                            ui.label("τ (fs)");
                            ui.add(egui::DragValue::new(&mut s.tau_fs).speed(1.0));
                        });
                    }
                    ui.label(
                        egui::RichText::new(
                            "Forces are numerical (finite-difference qchem energies), so keep \
                             systems small — runs on a background thread and is cost-capped.",
                        )
                        .weak()
                        .small(),
                    );

                    ui.add_space(6.0);
                    ui.add_enabled_ui(!running, |ui| {
                        if ui
                            .button(egui::RichText::new("▶ Run dynamics").strong())
                            .clicked()
                        {
                            start_run(s);
                        }
                    });

                    if running {
                        let done = s.run.as_ref().map(|r| r.progress.load(Ordering::Relaxed)).unwrap_or(0);
                        let total = s.run.as_ref().map(|r| r.total).unwrap_or(1).max(1);
                        ui.add(
                            egui::ProgressBar::new(done as f32 / total as f32)
                                .text(format!("step {done}/{total}")),
                        );
                        ui.ctx().request_repaint();
                    }
                    if !s.status.is_empty() {
                        ui.label(egui::RichText::new(&s.status).weak().small());
                    }
                    if let Some(e) = &s.error {
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    let ctx = ui.ctx().clone();
                    draw_results_and_playback(s, ui, &ctx, &mut to_push);
                });
        },
    );
    if close {
        app.show_reactdyn_workbench = false;
    }

    // The panel closure has released its &mut app borrow — now push the
    // selected frame's molecule into the shared 3-D viewport.
    if let Some((view, label)) = to_push {
        let mesh = molecule_view::ball_and_stick(&view, 0.28, 0.18);
        let _ = molecule_view::show_molecule(app, mesh, &label);
    }
}

/// Render the result summary + 3-D playback controls + energy plot, and
/// (when the displayed frame changes) stage the frame's molecule in
/// `to_push` for the caller to send to the viewport.
fn draw_results_and_playback(
    s: &mut ReactdynWorkbenchState,
    ui: &mut egui::Ui,
    ctx: &egui::Context,
    to_push: &mut Option<(ViewMolecule, String)>,
) {
    // Snapshot the energy series + frame count, dropping the s.last borrow
    // so the playback UI can mutate s.frame_idx / s.playing freely.
    let info = s.last.as_ref().map(|t| {
        let pe: Vec<[f64; 2]> = t
            .frames
            .iter()
            .map(|f| [f.time_fs, f.potential_hartree])
            .collect();
        let ke: Vec<[f64; 2]> = t
            .frames
            .iter()
            .map(|f| [f.time_fs, f.kinetic_hartree])
            .collect();
        let tot: Vec<[f64; 2]> = t
            .frames
            .iter()
            .map(|f| [f.time_fs, f.total_hartree()])
            .collect();
        (t.frames.len(), pe, ke, tot)
    });
    let Some((n, pe, ke, tot)) = info else {
        return;
    };
    if n == 0 {
        return;
    }

    ui.separator();
    ui.label(egui::RichText::new("Result").strong());
    if let Some(t) = &s.last {
        ui.label(egui::RichText::new(summarize(t)).monospace().small());
    }

    ui.add_space(4.0);
    ui.label(egui::RichText::new("Playback (3-D viewport)").strong());
    if s.frame_idx >= n {
        s.frame_idx = n - 1;
    }
    ui.horizontal(|ui| {
        if ui
            .button(if s.playing { "⏸ pause" } else { "▶ play" })
            .clicked()
        {
            s.playing = !s.playing;
        }
        if ui.button("⏮").on_hover_text("first frame").clicked() {
            s.frame_idx = 0;
            s.playing = false;
        }
        if ui.button("⏭").on_hover_text("last frame").clicked() {
            s.frame_idx = n - 1;
            s.playing = false;
        }
    });
    ui.add(egui::Slider::new(&mut s.frame_idx, 0..=n - 1).text("frame"));

    Plot::new("reactdyn_energy")
        .height(150.0)
        .legend(Legend::default())
        .show(ui, |pui| {
            pui.line(Line::new(PlotPoints::from(pe)).name("potential"));
            pui.line(Line::new(PlotPoints::from(ke)).name("kinetic"));
            pui.line(Line::new(PlotPoints::from(tot)).name("total"));
        });

    // Current-frame readout + stage the viewport mesh when the frame
    // changes (recompute bonds each frame → bonds visibly form/break).
    let cur = s.frame_idx;
    let need_push = s.last_pushed != Some(cur);
    if let Some(t) = &s.last {
        let f = &t.frames[cur];
        ui.label(
            egui::RichText::new(format!(
                "t = {:.2} fs    E_total = {:.6} Ha",
                f.time_fs,
                f.total_hartree()
            ))
            .monospace()
            .small(),
        );
        if need_push {
            *to_push = Some((
                view_molecule_for_frame(t, cur),
                format!("AIMD frame {cur}/{}", n - 1),
            ));
        }
    }
    if need_push {
        s.last_pushed = Some(cur);
    }
    if s.playing {
        s.frame_idx = (cur + 1) % n;
        ctx.request_repaint();
    }
}

/// Build a [`ViewMolecule`] for one trajectory frame: atoms at the frame's
/// positions (ångström) with bonds re-detected from the current geometry,
/// so a bond appears/disappears as it forms/breaks.
fn view_molecule_for_frame(traj: &Trajectory, idx: usize) -> ViewMolecule {
    let frame = &traj.frames[idx];
    let pos = frame.pos_angstrom();
    let atoms: Vec<ViewAtom> = pos
        .iter()
        .enumerate()
        .map(|(i, p)| {
            ViewAtom::new(
                [p[0] as f32, p[1] as f32, p[2] as f32],
                traj.system.elements[i].symbol(),
            )
        })
        .collect();
    let bonds = molecule_view::detect_bonds(&atoms);
    ViewMolecule { atoms, bonds }
}

/// A text summary of a finished trajectory — frame count, duration, and
/// the energy-conservation diagnostic (the honest correctness readout).
fn summarize(traj: &Trajectory) -> String {
    let frames = &traj.frames;
    let (e0, ef, t_end) = match (frames.first(), frames.last()) {
        (Some(a), Some(b)) => (a.total_hartree(), b.total_hartree(), b.time_fs),
        _ => (0.0, 0.0, 0.0),
    };
    let max_drift = frames
        .iter()
        .map(|f| (f.total_hartree() - e0).abs())
        .fold(0.0_f64, f64::max);
    format!(
        "atoms {} · {} frames · {:.2} fs\n\
         E(total): {:.6} → {:.6} Ha   (max drift {:.2e})",
        traj.system.n_atoms(),
        frames.len(),
        t_end,
        e0,
        ef,
        max_drift,
    )
}

/// The agent-bridge product for the reaction-dynamics workbench
/// (`show_3d{kind="reactdyn"}`).
///
/// Exposes **one representative frame** of the canonical reaction — the initial
/// geometry (frame 0) of the default **Water** preset, parsed straight from its
/// built-in XYZ via [`MolecularGeometry::from_xyz_str`] into a [`ViewMolecule`].
/// This is the same geometry an AIMD run starts from, but it needs **no** solve:
/// the numerical-gradient Born-Oppenheimer dynamics is deliberately *not* run
/// here, so the builder stays pure, deterministic and cheap (no background
/// thread, no qchem cost). The frame is meshed as a colour-aware ball-and-stick
/// (CPK by element, reusing the molecule view's
/// [`molecule_view::ball_and_stick_colored`]) and promoted to a `Tri3`
/// [`valenx_mesh::Mesh`] the tile renders coloured. The readout names the
/// reaction system and its atom count.
pub(crate) fn reactdyn_product() -> crate::WorkspaceProduct {
    // The representative frame: the Water preset's initial geometry (frame 0).
    let mol = (|| -> Option<ViewMolecule> {
        let xyz = Preset::Water.xyz()?;
        let geom = MolecularGeometry::from_xyz_str(xyz).ok()?;
        if geom.atoms.is_empty() {
            return None;
        }
        let atoms: Vec<ViewAtom> = geom
            .atoms
            .iter()
            .map(|a| {
                let p = a.position_angstrom();
                ViewAtom::new([p[0] as f32, p[1] as f32, p[2] as f32], a.element.symbol())
            })
            .collect();
        let bonds = molecule_view::detect_bonds(&atoms);
        Some(ViewMolecule { atoms, bonds })
    })()
    .unwrap_or_default();

    let (soup, per_tri_colors) = molecule_view::ball_and_stick_colored(&mol, 0.28, 0.18);
    let mesh = crate::products_registry::mesh_from_triangle_soup(&soup, "valenx-reactdyn-frame");
    let vertex_colors = crate::products_registry::per_triangle_to_vertex_colors(&per_tri_colors);
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<reactdyn>/frame-0");
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    let lines = vec![
        "reaction dynamics: H₂O (Born-Oppenheimer AIMD)".to_string(),
        format!("representative frame 0 · {} atoms", mol.atoms.len()),
        "initial geometry — run the workbench for the full trajectory".to_string(),
    ];
    crate::WorkspaceProduct {
        title: "Reaction dynamics (frame)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: Some(vertex_colors),
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
    fn default_h2_inputs_run_end_to_end() {
        let mut s = ReactdynWorkbenchState {
            n_steps: 10, // small for a fast headless test
            ..Default::default()
        };
        let (system, controls) = build_aimd_inputs(&s).expect("inputs should parse");
        assert_eq!(system.n_atoms(), 2);
        let traj = AimdEngine
            .run(&system, &controls, &mut |_| {})
            .expect("AIMD run should succeed");
        assert_eq!(traj.frames.len(), controls.n_steps + 1);
        s.last = Some(traj);
        assert!(summarize(s.last.as_ref().unwrap()).contains("max drift"));
    }

    #[test]
    fn bad_custom_xyz_fails_loud() {
        let s = ReactdynWorkbenchState {
            preset: Preset::Custom,
            xyz: "not a geometry".into(),
            ..Default::default()
        };
        assert!(build_aimd_inputs(&s).is_err());
    }

    #[test]
    fn frame_view_has_atoms_and_detects_the_h2_bond() {
        let s = ReactdynWorkbenchState {
            n_steps: 3,
            ..Default::default()
        };
        let (system, controls) = build_aimd_inputs(&s).unwrap();
        let traj = AimdEngine.run(&system, &controls, &mut |_| {}).unwrap();
        let view = view_molecule_for_frame(&traj, 0);
        assert_eq!(view.atoms.len(), 2);
        // H2 at ~0.95 Å is bonded → exactly one bond.
        assert_eq!(
            view.bonds.len(),
            1,
            "expected one H-H bond, got {:?}",
            view.bonds
        );
    }

    #[test]
    fn qmmm_inputs_build_and_run_with_solvent_shell() {
        let s = ReactdynWorkbenchState {
            backend: Backend::QmMm,
            n_solvent: 8,
            n_steps: 4,
            ..Default::default()
        };
        let (sys, controls) = build_qmmm_inputs(&s).expect("qmmm inputs build");
        assert_eq!(sys.n_qm(), 2); // H2 default
        assert_eq!(sys.n_mm(), 8); // the solvent shell
        let traj = QmMmEngine
            .run(&sys, &controls, &mut |_| {})
            .expect("qmmm run");
        assert_eq!(traj.system.n_atoms(), 10); // 2 QM + 8 MM atoms in 3-D
    }

    #[test]
    fn electrostatic_mode_runs_via_the_workbench() {
        let s = ReactdynWorkbenchState {
            backend: Backend::QmMm,
            embedding: Embedding::Electrostatic,
            n_solvent: 4,
            n_steps: 3,
            ..Default::default()
        };
        let (sys, controls) = build_qmmm_inputs(&s).expect("electrostatic inputs");
        let traj = QmMmEngine
            .run(&sys, &controls, &mut |_| {})
            .expect("electrostatic run");
        assert_eq!(traj.system.n_atoms(), 6); // 2 QM + 4 MM
    }

    #[test]
    fn reactive_cluster_runs_via_the_workbench() {
        let s = ReactdynWorkbenchState {
            backend: Backend::Reactive,
            n_cluster: 6,
            n_steps: 5,
            dt_fs: 0.2,
            ..Default::default()
        };
        let (sys, controls) = build_reactive_inputs(&s).expect("reactive inputs");
        assert_eq!(sys.n_atoms(), 6);
        let traj = ReactiveEngine
            .run(&sys, &controls, &mut |_| {})
            .expect("reactive run");
        assert_eq!(traj.system.n_atoms(), 6);
        assert_eq!(traj.frames.len(), 6);
    }
}
