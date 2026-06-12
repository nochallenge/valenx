//! Panel 9 — **Quantum Chemistry** (`valenx-qchem`).
//!
//! Enter a molecular geometry (XYZ), run a restricted or unrestricted
//! Hartree-Fock SCF in a built-in Gaussian basis, and read off the
//! total energy, molecular-orbital energies, the HOMO-LUMO gap and
//! Mulliken atomic charges — all native `valenx-qchem` calls.

use eframe::egui;

use valenx_qchem::dft::{Functional, GridQuality};
use valenx_qchem::driver::{run_dft, run_mp2, run_rhf, run_uhf};
use valenx_qchem::geometry::MolecularGeometry;
use valenx_qchem::scf::rhf::ScfSettings;

use super::common;
use crate::ValenxApp;

/// Electronic-structure method choice.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Method {
    #[default]
    Rhf,
    Uhf,
    /// Restricted Kohn-Sham DFT (functional selected via [`Xc`]).
    Dft,
    /// RHF reference + MP2 correlation correction.
    Mp2,
}

/// Exchange-correlation functional for DFT — maps to
/// [`valenx_qchem::dft::Functional`].
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum Xc {
    Lda,
    Pbe,
    #[default]
    B3lyp,
}

impl Xc {
    fn to_functional(self) -> Functional {
        match self {
            Xc::Lda => Functional::Lda,
            Xc::Pbe => Functional::Pbe,
            Xc::B3lyp => Functional::B3lyp,
        }
    }
}

/// Snapshot of every editable input the Qchem panel owns.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub(crate) struct QchemSnapshot {
    pub(crate) method: Method,
    pub(crate) functional: Xc,
    pub(crate) basis: String,
    pub(crate) geometry: String,
}

/// Form + result state for the Quantum Chemistry panel.
pub struct QchemPanel {
    method: Method,
    /// Exchange-correlation functional (used when `method == Dft`).
    functional: Xc,
    /// Basis-set name (one of the built-in sets).
    basis: String,
    /// XYZ geometry text (the comment line may carry `charge mult`).
    geometry: String,
    error: Option<String>,
    result: String,
    /// Undo / redo over the geometry text + method + basis.
    history: crate::undo::History<QchemSnapshot>,
}

impl QchemPanel {
    fn snapshot(&self) -> QchemSnapshot {
        QchemSnapshot {
            method: self.method,
            functional: self.functional,
            basis: self.basis.clone(),
            geometry: self.geometry.clone(),
        }
    }
    fn restore(&mut self, s: QchemSnapshot) {
        self.method = s.method;
        self.functional = s.functional;
        self.basis = s.basis;
        self.geometry = s.geometry;
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

impl Default for QchemPanel {
    fn default() -> Self {
        QchemPanel {
            method: Method::Rhf,
            functional: Xc::default(),
            basis: "STO-3G".to_string(),
            // Water at a reasonable geometry — closed-shell, RHF-ready.
            geometry: "3\nwater  charge 0  mult 1\n\
                       O  0.000000  0.000000  0.117300\n\
                       H  0.000000  0.757200 -0.469200\n\
                       H  0.000000 -0.757200 -0.469200\n"
                .to_string(),
            error: None,
            result: String::new(),
            history: crate::undo::History::new(),
        }
    }
}

/// Render the Quantum Chemistry panel.
pub fn draw(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let p = &mut app.genetics.qchem;

    common::section(ui, "Molecular geometry (XYZ)");
    ui.label(
        egui::RichText::new("comment line may carry `charge N mult M`")
            .weak()
            .small(),
    );
    ui.horizontal(|ui| {
        if ui.small_button("Load XYZ…").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("XYZ", &["xyz", "txt"])
                .pick_file()
            {
                // Round-21 H1: see biostruct loader.
                match valenx_core::io_caps::read_capped_to_string(
                    &path,
                    valenx_core::io_caps::MAX_GENETICS_FILE_BYTES as usize,
                ) {
                    Ok(t) => p.geometry = t,
                    Err(e) => p.error = Some(format!("read: {e}")),
                }
            }
        }
    });
    ui.add(
        egui::TextEdit::multiline(&mut p.geometry)
            .id_source("qchem_geometry")
            .font(egui::TextStyle::Monospace)
            .desired_width(f32::INFINITY)
            .desired_rows(6),
    );
    ui.separator();

    common::section(ui, "Method + basis");
    ui.horizontal_wrapped(|ui| {
        ui.radio_value(&mut p.method, Method::Rhf, "RHF")
            .on_hover_text("Restricted Hartree-Fock (closed-shell)");
        ui.radio_value(&mut p.method, Method::Uhf, "UHF")
            .on_hover_text("Unrestricted Hartree-Fock (open-shell)");
        ui.radio_value(&mut p.method, Method::Dft, "DFT")
            .on_hover_text("Restricted Kohn-Sham density-functional theory");
        ui.radio_value(&mut p.method, Method::Mp2, "MP2")
            .on_hover_text("RHF reference + MP2 correlation (closed-shell)");
    });
    if p.method == Method::Dft {
        ui.horizontal(|ui| {
            ui.label("Functional:");
            ui.radio_value(&mut p.functional, Xc::Lda, "LDA");
            ui.radio_value(&mut p.functional, Xc::Pbe, "PBE");
            ui.radio_value(&mut p.functional, Xc::B3lyp, "B3LYP");
        });
    }
    ui.horizontal(|ui| {
        ui.label("Basis set:");
        egui::ComboBox::from_id_source("qchem_basis")
            .selected_text(&p.basis)
            .show_ui(ui, |ui| {
                for b in ["STO-3G", "3-21G", "6-31G", "6-31G*"] {
                    ui.selectable_value(&mut p.basis, b.to_string(), b);
                }
            });
    });

    if common::run_button(ui, "Run SCF") {
        let snap = p.snapshot();
        p.history.record(snap);
        run_scf(p);
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
    if !p.result.is_empty() {
        ui.separator();
        common::section(ui, "Result");
        common::mono_output(ui, "qchem_result", &p.result, 18);
    }
}

/// Run the Hartree-Fock SCF — extracted from the button closure so it
/// is callable from the headless UI tests.
fn run_scf(p: &mut QchemPanel) {
    p.error = None;
    match MolecularGeometry::from_xyz_str(&p.geometry) {
        Ok(geom) => {
            let settings = ScfSettings::default();
            let report = match p.method {
                Method::Rhf => run_rhf(&geom, &p.basis, settings),
                Method::Uhf => run_uhf(&geom, &p.basis, settings),
                Method::Dft => run_dft(
                    &geom,
                    &p.basis,
                    p.functional.to_functional(),
                    GridQuality::default(),
                    settings,
                ),
                Method::Mp2 => run_mp2(&geom, &p.basis, settings),
            };
            match report {
                Ok(r) => {
                    let method_label = match &r.dft {
                        Some(d) => format!("DFT / {}", d.functional.label()),
                        None => r.method.label().to_string(),
                    };
                    let mut out = format!(
                        "method         : {}\nbasis          : {}\n\
                             SCF iterations : {}\n\n\
                             total energy   : {:.8} Ha\n",
                        method_label, r.basis_name, r.scf_iterations, r.total_energy,
                    );
                    // MP2: show the HF reference + correlation correction.
                    if let Some(corr) = r.correlation_energy {
                        out.push_str(&format!(
                            "Hartree-Fock E : {:.8} Ha\n\
                                 correlation E  : {:.8} Ha\n",
                            r.hartree_fock_energy, corr,
                        ));
                    }
                    // DFT: show the exchange-correlation energy + a grid
                    // sanity check (integrated electron count).
                    if let Some(d) = &r.dft {
                        out.push_str(&format!(
                            "XC energy Exc  : {:.8} Ha\n\
                                 grid electrons : {:.4}  (≈ N electrons)\n",
                            d.xc_energy, d.grid_electron_count,
                        ));
                    }
                    out.push_str(&format!(
                        "nuclear rep.   : {:.8} Ha\n\
                             dipole moment  : {:.4} D\n",
                        r.nuclear_repulsion,
                        r.dipole.magnitude_debye(),
                    ));
                    match r.homo_lumo_gap {
                        Some(g) => out.push_str(&format!(
                            "HOMO-LUMO gap  : {:.4} Ha  ({:.2} eV)\n",
                            g,
                            g * 27.211_386,
                        )),
                        None => out.push_str("HOMO-LUMO gap  : (n/a)\n"),
                    }
                    if let Some(s2) = r.s_squared {
                        out.push_str(&format!("⟨S²⟩           : {s2:.4}\n"));
                    }

                    out.push_str("\n-- molecular orbitals --\n");
                    let orbs = &r.orbitals.orbitals;
                    let homo = r.orbitals.homo_index;
                    for (i, mo) in orbs.iter().enumerate() {
                        let tag = match homo {
                            Some(h) if i == h => "  <- HOMO",
                            Some(h) if i == h + 1 => "  <- LUMO",
                            _ => "",
                        };
                        out.push_str(&format!(
                            "  MO {:<3} {:>10.5} Ha  occ {:.1}{}\n",
                            i + 1,
                            mo.energy,
                            mo.occupation,
                            tag,
                        ));
                    }

                    out.push_str("\n-- Mulliken atomic charges --\n");
                    for (i, q) in r.partial_charges.iter().enumerate() {
                        out.push_str(&format!("  atom {:<3} {:>8.4} e\n", i + 1, q));
                    }
                    p.result = out;
                }
                Err(e) => p.error = Some(e.to_string()),
            }
        }
        Err(e) => p.error = Some(e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_geometry_parses() {
        let p = QchemPanel::default();
        let geom = MolecularGeometry::from_xyz_str(&p.geometry).expect("water XYZ must parse");
        assert_eq!(geom.atoms.len(), 3);
    }

    #[test]
    fn water_rhf_converges() {
        // STO-3G RHF on water is small and fast — a real end-to-end
        // SCF check that the panel's wiring is sound.
        let p = QchemPanel::default();
        let geom = MolecularGeometry::from_xyz_str(&p.geometry).unwrap();
        let r = run_rhf(&geom, "STO-3G", ScfSettings::default()).unwrap();
        assert!(r.total_energy < 0.0, "HF energy must be negative");
    }

    #[test]
    fn water_dft_b3lyp_converges_and_reports_xc() {
        // The native Kohn-Sham DFT path the panel now surfaces: a real
        // end-to-end B3LYP SCF on water, with the DFT-specific report.
        let p = QchemPanel::default();
        let geom = MolecularGeometry::from_xyz_str(&p.geometry).unwrap();
        let r = run_dft(
            &geom,
            "STO-3G",
            Functional::B3lyp,
            GridQuality::default(),
            ScfSettings::default(),
        )
        .unwrap();
        assert!(r.total_energy < 0.0, "KS energy must be negative");
        assert!(r.dft.is_some(), "DFT report must carry DftInfo");
    }
}

/// Headless egui UI-logic tests for the Quantum Chemistry panel.
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
        app.genetics.active = GeneticsPanel::QuantumChemistry;
        app
    }

    #[test]
    fn draws_both_methods_without_panic() {
        for method in [Method::Rhf, Method::Uhf] {
            let mut app = app_with_panel();
            app.genetics.qchem.method = method;
            draw_headless(&mut app);
        }
    }

    #[test]
    fn draws_post_run_and_error_states_without_panic() {
        let mut app = app_with_panel();
        app.genetics.qchem.result = "method : Rhf\ntotal energy : -74.96 Ha\n".to_string();
        draw_headless(&mut app);
        let mut app = app_with_panel();
        app.genetics.qchem.error = Some("geometry parse failed".to_string());
        draw_headless(&mut app);
    }

    #[test]
    fn run_scf_converges_water_rhf() {
        // STO-3G RHF on the default water geometry → the real
        // valenx-qchem SCF driver produces a correctly-formatted
        // result with a negative total energy.
        let mut p = QchemPanel::default();
        run_scf(&mut p);
        assert!(p.error.is_none(), "SCF errored: {:?}", p.error);
        assert!(p.result.contains("total energy"));
        assert!(p.result.contains("molecular orbitals"));
        assert!(p.result.contains("Mulliken"));
    }

    #[test]
    fn run_scf_surfaces_error_on_malformed_geometry() {
        // Junk text is not a valid XYZ geometry — the panel must
        // surface an error rather than panicking.
        let mut p = QchemPanel {
            geometry: "not an xyz file at all".to_string(),
            ..QchemPanel::default()
        };
        run_scf(&mut p);
        assert!(p.error.is_some(), "SCF should error on malformed geometry");
        assert!(p.result.is_empty());
    }
}
