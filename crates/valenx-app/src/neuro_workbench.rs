//! Neural-interface (BCI stimulation) workbench.
//!
//! A right-side panel driving `valenx-neuro`: place a stimulating electrode in
//! tissue with a bundle of axons, set the current, and Run — the panel shows
//! which axons are recruited, the recruitment curve, the electrode-impedance
//! Bode plot, the tissue temperature rise, and a cross-section schematic.
//! Compute is synchronous (a few seconds); a background runner is future work.

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use valenx_neuro::{
    recruitment_curve, solve_point_heat, stimulate, Axon, Cpe, ElectrodeImpedance, NeuroError,
    Scene, TissueGrid,
};

use crate::ValenxApp;

/// Persistent state for the neural-interface workbench.
pub struct NeuroWorkbenchState {
    electrode_ua: f64,
    electrode_radius_um: f64,
    sigma_s_m: f64,
    n_axons: usize,
    depth_mm: f64,
    spread_mm: f64,
    results: Option<NeuroResults>,
    error: Option<String>,
}

impl Default for NeuroWorkbenchState {
    fn default() -> Self {
        Self {
            electrode_ua: 300.0,
            electrode_radius_um: 50.0,
            sigma_s_m: 0.2,
            n_axons: 8,
            depth_mm: 0.5,
            spread_mm: 2.0,
            results: None,
            error: None,
        }
    }
}

struct NeuroResults {
    recruited_fraction: f64,
    fired: Vec<bool>,
    depths_mm: Vec<f64>,
    recruitment_curve: Vec<(f64, f64)>,
    impedance_bode: Vec<(f64, f64)>,
    access_resistance_ohm: f64,
    dt_at_1mm_k: f64,
}

/// Run one stimulation + a recruitment sweep for the current settings.
fn run_neuro(s: &NeuroWorkbenchState) -> Result<NeuroResults, NeuroError> {
    if s.n_axons == 0 {
        return Err(NeuroError::Invalid {
            reason: "add at least one axon".to_string(),
        });
    }
    let sigma = s.sigma_s_m.max(1.0e-3);
    let grid = TissueGrid::cube(40.0, 21, sigma);
    let axons: Vec<Axon> = (0..s.n_axons)
        .map(|i| {
            let frac = if s.n_axons > 1 {
                i as f64 / (s.n_axons as f64 - 1.0)
            } else {
                0.0
            };
            Axon::squid_at(s.depth_mm.max(0.1) + frac * s.spread_mm.max(0.0))
        })
        .collect();
    let depths_mm: Vec<f64> = axons.iter().map(|a| a.depth_mm).collect();
    let scene = Scene { grid, axons };

    let mag = s.electrode_ua.abs();
    let rec = stimulate(&scene, -mag)?;
    let fired = rec.fired().to_vec();
    let recruited_fraction = rec.recruited_fraction();

    let curve = recruitment_curve(&scene, &[10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0])?;

    let imp = ElectrodeImpedance::disk(s.electrode_radius_um.max(1.0), sigma, Cpe::default());
    let access_resistance_ohm = imp.access_resistance_ohm();
    let i_amp = mag * 1.0e-6;
    let power_w = i_amp * i_amp * access_resistance_ohm; // I²R Joule heating
    let bio = solve_point_heat(&scene.grid, 0.5, power_w)?;
    let dt_at_1mm_k = bio.delta_t_k_at_radius_x(1.0);

    let impedance_bode: Vec<(f64, f64)> = [1.0, 10.0, 100.0, 1.0e3, 1.0e4, 1.0e5]
        .iter()
        .map(|&f| (f, imp.magnitude_ohm(f)))
        .collect();

    Ok(NeuroResults {
        recruited_fraction,
        fired,
        depths_mm,
        recruitment_curve: curve,
        impedance_bode,
        access_resistance_ohm,
        dt_at_1mm_k,
    })
}

/// Draw the neural-interface workbench (a no-op unless toggled on via
/// View → Neural Interface).
pub fn draw_neuro_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_neuro_workbench {
        return;
    }
    egui::SidePanel::right("valenx_neuro_workbench")
        .resizable(true)
        .default_width(380.0)
        .width_range(320.0..=640.0)
        .show(ctx, |ui| {
            ui.heading("Neural Interface");
            ui.label(
                egui::RichText::new("native BCI stimulation · valenx-neuro")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.neuro;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Electrode & tissue").strong());
                    ui.horizontal(|ui| {
                        ui.label("current (µA, cathodic)");
                        ui.add(egui::DragValue::new(&mut s.electrode_ua).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("electrode radius (µm)");
                        ui.add(egui::DragValue::new(&mut s.electrode_radius_um).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("tissue σ (S/m)");
                        ui.add(egui::DragValue::new(&mut s.sigma_s_m).speed(0.01));
                    });
                    ui.label(egui::RichText::new("Axon bundle").strong());
                    ui.horizontal(|ui| {
                        ui.label("axons");
                        ui.add(egui::DragValue::new(&mut s.n_axons).speed(0.3));
                    });
                    ui.horizontal(|ui| {
                        ui.label("nearest depth (mm)");
                        ui.add(egui::DragValue::new(&mut s.depth_mm).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("spread (mm)");
                        ui.add(egui::DragValue::new(&mut s.spread_mm).speed(0.1));
                    });
                    ui.separator();
                    if ui.button("▶ Run stimulation").clicked() {
                        match run_neuro(s) {
                            Ok(r) => {
                                s.results = Some(r);
                                s.error = None;
                            }
                            Err(e) => {
                                s.error = Some(e.to_string());
                                s.results = None;
                            }
                        }
                    }
                    if let Some(e) = &s.error {
                        ui.colored_label(egui::Color32::RED, e);
                    }
                    if let Some(r) = &s.results {
                        ui.separator();
                        let n_fired = r.fired.iter().filter(|&&f| f).count();
                        ui.label(
                            egui::RichText::new(format!(
                                "recruited {:.0}%  ({n_fired}/{})\naccess resistance {:.1} kΩ\ntissue ΔT @1 mm {:.3} K",
                                r.recruited_fraction * 100.0,
                                r.fired.len(),
                                r.access_resistance_ohm / 1000.0,
                                r.dt_at_1mm_k,
                            ))
                            .monospace()
                            .small(),
                        );
                        ui.add_space(4.0);
                        ui.label(egui::RichText::new("Recruitment vs current (µA)").strong());
                        Plot::new("neuro_recruitment")
                            .height(130.0)
                            .legend(Legend::default())
                            .show(ui, |pui| {
                                let pts: Vec<[f64; 2]> =
                                    r.recruitment_curve.iter().map(|&(c, f)| [c, f]).collect();
                                pui.line(Line::new(PlotPoints::from(pts)).name("fraction"));
                            });
                        ui.label(egui::RichText::new("Electrode |Z| Bode (log–log)").strong());
                        Plot::new("neuro_impedance").height(130.0).show(ui, |pui| {
                            let pts: Vec<[f64; 2]> = r
                                .impedance_bode
                                .iter()
                                .map(|&(f, z)| [f.log10(), z.log10()])
                                .collect();
                            pui.line(Line::new(PlotPoints::from(pts)).name("|Z|"));
                        });
                        ui.label(
                            egui::RichText::new("Cross-section (● electrode, — axons)")
                                .weak()
                                .small(),
                        );
                        draw_schematic(ui, r);
                    }
                });
        });
}

/// A 2-D cross-section: the electrode at top, axons at their depths, coloured
/// green when recruited.
fn draw_schematic(ui: &mut egui::Ui, r: &NeuroResults) {
    let (resp, painter) =
        ui.allocate_painter(egui::vec2(ui.available_width(), 120.0), egui::Sense::hover());
    let rect = resp.rect;
    let cx = rect.center().x;
    let top = rect.top() + 12.0;
    painter.circle_filled(egui::pos2(cx, top), 5.0, egui::Color32::YELLOW);
    let max_depth = r.depths_mm.iter().cloned().fold(0.5_f64, f64::max) + 0.5;
    for (&d, &fired) in r.depths_mm.iter().zip(&r.fired) {
        let y = top + 14.0 + (d / max_depth) as f32 * (rect.height() - 28.0);
        let color = if fired {
            egui::Color32::from_rgb(80, 220, 120)
        } else {
            egui::Color32::GRAY
        };
        painter.line_segment(
            [
                egui::pos2(rect.left() + 18.0, y),
                egui::pos2(rect.right() - 18.0, y),
            ],
            egui::Stroke::new(2.0, color),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_neuro_default_produces_sane_results() {
        let s = NeuroWorkbenchState::default();
        let r = run_neuro(&s).expect("run");
        assert!((0.0..=1.0).contains(&r.recruited_fraction));
        assert_eq!(r.fired.len(), 8);
        assert!(r.access_resistance_ohm > 0.0);
        assert!(r.dt_at_1mm_k >= 0.0);
        let fr: Vec<f64> = r.recruitment_curve.iter().map(|&(_, f)| f).collect();
        for w in fr.windows(2) {
            assert!(w[1] >= w[0], "recruitment curve must not decrease: {fr:?}");
        }
    }

    #[test]
    fn run_neuro_zero_axons_errors() {
        let s = NeuroWorkbenchState {
            n_axons: 0,
            ..Default::default()
        };
        assert!(run_neuro(&s).is_err());
    }
}
