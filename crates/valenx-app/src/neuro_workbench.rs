//! Neural-interface (BCI stimulation) workbench.
//!
//! A right-side panel driving `valenx-neuro`: place a stimulating electrode in
//! tissue with a bundle of axons, set the current, and Run — the panel shows
//! which axons are recruited, the recruitment curve, the electrode-impedance
//! Bode plot, the tissue temperature rise, the recorded extracellular spike
//! (EAP), and a cross-section schematic.
//! Compute is synchronous (a few seconds); a background runner is future work.

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};
use nalgebra::Vector3;

use valenx_neuro::{
    activation_radius, charge_density, is_safe, max_safe_charge_per_phase, recruitment_curve,
    shannon_k, solve_point_heat, stimulate, threshold_current, Axon, Cpe, ElectrodeImpedance,
    ExtracellularRecorder, NeuroError, Scene, TissueGrid,
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
    /// Stimulus pulse width (µs) — charge per phase = current × width.
    pulse_width_us: f64,
    /// Shannon damage-line limit (k-value); points below it are safe.
    k_limit: f64,
    /// Current–distance: at-electrode threshold I₀ (µA).
    cd_i0_ua: f64,
    /// Current–distance constant k (µA/µm²): `I_th = I₀ + k·r²`.
    cd_k_ua_per_um2: f64,
    /// Cable theory: fiber diameter (µm).
    fiber_diameter_um: f64,
    /// Cable theory: specific membrane resistance R_m (Ω·cm²).
    r_m_ohm_cm2: f64,
    /// Cable theory: axial resistivity R_i (Ω·cm).
    r_i_ohm_cm: f64,
    /// Cable theory: specific membrane capacitance C_m (µF/cm²).
    c_m_uf_cm2: f64,
    /// Conduction: Hursh factor for myelinated velocity (m/s per µm).
    k_myel_m_per_s_per_um: f64,
    /// Conduction: cable scaling for unmyelinated velocity (m/s per √µm).
    k_unmyel_m_per_s_per_sqrt_um: f64,
    /// Conduction: propagation distance for the latency readout (m).
    conduction_distance_m: f64,
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
            pulse_width_us: 200.0,
            k_limit: 1.85,
            cd_i0_ua: 2.0,
            cd_k_ua_per_um2: 0.001,
            fiber_diameter_um: 1.0,
            r_m_ohm_cm2: 10000.0,
            r_i_ohm_cm: 100.0,
            c_m_uf_cm2: 1.0,
            k_myel_m_per_s_per_um: valenx_neuro::HURSH_FACTOR_M_PER_S_PER_UM,
            k_unmyel_m_per_s_per_sqrt_um: 2.0,
            conduction_distance_m: 0.05,
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
    eap_uv: Vec<f64>,
    eap_dt_ms: f64,
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

    // Forward recording: the extracellular spike (EAP) a nearby electrode would
    // see from a representative axon at the bundle's nearest depth.
    let recorder = ExtracellularRecorder::new(sigma, 100.0, 238.0, 35.4);
    let eap = recorder.record(
        200,
        Vector3::new(10.0e-3, s.depth_mm.max(0.1) * 1.0e-3, 0.0),
    );

    Ok(NeuroResults {
        recruited_fraction,
        fired,
        depths_mm,
        recruitment_curve: curve,
        impedance_bode,
        access_resistance_ohm,
        dt_at_1mm_k,
        eap_uv: eap.eap_uv,
        eap_dt_ms: eap.dt_ms,
    })
}

/// Charge-injection safety readout for the current electrode + pulse, via the
/// `valenx-neuro` Shannon model.
struct SafetyReadout {
    q_uc: f64,
    area_cm2: f64,
    density_uc_cm2: f64,
    k_value: f64,
    safe: bool,
    max_safe_q_uc: f64,
}

/// Compute the charge-injection safety of the current settings. Charge per
/// phase `Q = I·PW` (µA·µs → µC); electrode area from the disk radius.
fn charge_safety(s: &NeuroWorkbenchState) -> SafetyReadout {
    // µA × µs = 1e-12 C = 1e-6 µC.
    let q_uc = s.electrode_ua.abs() * s.pulse_width_us.max(0.0) * 1.0e-6;
    let r_cm = s.electrode_radius_um.max(1.0) * 1.0e-4; // µm → cm
    let area_cm2 = std::f64::consts::PI * r_cm * r_cm;
    SafetyReadout {
        q_uc,
        area_cm2,
        density_uc_cm2: charge_density(q_uc, area_cm2),
        k_value: shannon_k(q_uc, area_cm2),
        safe: is_safe(q_uc, area_cm2, s.k_limit),
        max_safe_q_uc: max_safe_charge_per_phase(area_cm2, s.k_limit),
    }
}

/// Current–distance selectivity readout for the current settings.
struct CdReadout {
    activation_radius_um: f64,
    threshold_at_depth_ua: f64,
    reaches_nearest: bool,
}

/// Apply the Stoney current–distance law to the current stimulus + geometry:
/// the activation radius for `|I|`, and the threshold to reach the nearest
/// axon at `depth_mm`.
fn current_distance_readout(s: &NeuroWorkbenchState) -> CdReadout {
    let i = s.electrode_ua.abs();
    let depth_um = s.depth_mm.max(0.0) * 1000.0;
    let threshold_at_depth_ua = threshold_current(s.cd_i0_ua, s.cd_k_ua_per_um2, depth_um);
    CdReadout {
        activation_radius_um: activation_radius(s.cd_i0_ua, s.cd_k_ua_per_um2, i),
        threshold_at_depth_ua,
        reaches_nearest: i >= threshold_at_depth_ua,
    }
}

/// Draw the neural-interface workbench (a no-op unless toggled on via
/// View → Neural Interface).
pub fn draw_neuro_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_neuro_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(app, ctx, "valenx_neuro_workbench", "Neural Interface", |app, ui| {
            ui.label(egui::RichText::new("native BCI stimulation · valenx-neuro").weak().small());
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
                    ui.label(egui::RichText::new("Charge-injection safety (Shannon)").strong());
                    ui.horizontal(|ui| {
                        ui.label("pulse width (µs)");
                        ui.add(egui::DragValue::new(&mut s.pulse_width_us).speed(5.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("Shannon k limit");
                        ui.add(egui::DragValue::new(&mut s.k_limit).speed(0.05));
                    });
                    {
                        let sf = charge_safety(s);
                        ui.label(
                            egui::RichText::new(format!(
                                "Q {:.3} µC/ph · area {:.2e} cm² · density {:.0} µC/cm² · k = {:.2}\nmax safe Q (k={:.2}): {:.3} µC/ph",
                                sf.q_uc,
                                sf.area_cm2,
                                sf.density_uc_cm2,
                                sf.k_value,
                                s.k_limit,
                                sf.max_safe_q_uc,
                            ))
                            .monospace()
                            .small(),
                        );
                        let (txt, col) = if sf.safe {
                            (
                                "● SAFE — below the Shannon limit",
                                egui::Color32::from_rgb(80, 220, 120),
                            )
                        } else {
                            (
                                "● UNSAFE — above the Shannon limit",
                                egui::Color32::from_rgb(220, 90, 90),
                            )
                        };
                        ui.colored_label(col, txt);
                    }

                    ui.separator();
                    ui.label(
                        egui::RichText::new("Current–distance selectivity (Stoney)").strong(),
                    );
                    ui.horizontal(|ui| {
                        ui.label("at-electrode threshold I₀ (µA)");
                        ui.add(egui::DragValue::new(&mut s.cd_i0_ua).speed(0.5));
                    });
                    ui.horizontal(|ui| {
                        ui.label("current–distance k (µA/µm²)");
                        ui.add(egui::DragValue::new(&mut s.cd_k_ua_per_um2).speed(0.0005));
                    });
                    {
                        let cd = current_distance_readout(s);
                        ui.label(
                            egui::RichText::new(format!(
                                "activation radius {:.0} µm ({:.2} mm) at I = {:.0} µA\nthreshold at nearest axon ({:.1} mm): {:.0} µA",
                                cd.activation_radius_um,
                                cd.activation_radius_um / 1000.0,
                                s.electrode_ua.abs(),
                                s.depth_mm,
                                cd.threshold_at_depth_ua,
                            ))
                            .monospace()
                            .small(),
                        );
                        let (txt, col) = if cd.reaches_nearest {
                            (
                                "● nearest axon within activation radius",
                                egui::Color32::from_rgb(80, 220, 120),
                            )
                        } else {
                            (
                                "● nearest axon beyond activation radius",
                                egui::Color32::from_rgb(220, 170, 80),
                            )
                        };
                        ui.colored_label(col, txt);
                    }

                    ui.separator();
                    ui.label(egui::RichText::new("Cable theory (passive fiber)").strong());
                    ui.horizontal(|ui| {
                        ui.label("fiber diameter (µm)");
                        ui.add(egui::DragValue::new(&mut s.fiber_diameter_um).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("R_m (Ω·cm²)");
                        ui.add(egui::DragValue::new(&mut s.r_m_ohm_cm2).speed(100.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("R_i (Ω·cm)");
                        ui.add(egui::DragValue::new(&mut s.r_i_ohm_cm).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("C_m (µF/cm²)");
                        ui.add(egui::DragValue::new(&mut s.c_m_uf_cm2).speed(0.05));
                    });
                    {
                        // µm → cm for the diameter; µF/cm² → F/cm² for the capacitance.
                        let d_cm = s.fiber_diameter_um * 1.0e-4;
                        let lambda_cm =
                            valenx_neuro::length_constant_cm(d_cm, s.r_m_ohm_cm2, s.r_i_ohm_cm);
                        let tau_s =
                            valenx_neuro::time_constant_s(s.r_m_ohm_cm2, s.c_m_uf_cm2 * 1.0e-6);
                        let f_c = valenx_neuro::membrane_cutoff_frequency(tau_s);
                        let r_inf = valenx_neuro::semi_infinite_input_resistance(
                            d_cm,
                            s.r_m_ohm_cm2,
                            s.r_i_ohm_cm,
                        );
                        ui.label(
                            egui::RichText::new(format!(
                                "space constant λ {:.4} cm ({:.0} µm)\nmembrane τ {:.2} ms · cutoff f_c {:.1} Hz · R_∞ {:.1} MΩ",
                                lambda_cm,
                                lambda_cm * 1.0e4,
                                tau_s * 1000.0,
                                f_c,
                                r_inf / 1.0e6,
                            ))
                            .monospace()
                            .small(),
                        );
                    }

                    ui.separator();
                    ui.label(egui::RichText::new("Conduction (diameter-driven)").strong());
                    ui.horizontal(|ui| {
                        ui.label("k_myelinated (m/s per µm)");
                        ui.add(egui::DragValue::new(&mut s.k_myel_m_per_s_per_um).speed(0.1));
                    });
                    ui.horizontal(|ui| {
                        ui.label("k_unmyelinated (m/s per √µm)");
                        ui.add(
                            egui::DragValue::new(&mut s.k_unmyel_m_per_s_per_sqrt_um).speed(0.05),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("distance for latency (m)");
                        ui.add(egui::DragValue::new(&mut s.conduction_distance_m).speed(0.005));
                    });
                    {
                        // fiber_diameter_um (cable-theory input above) is already in µm — exactly
                        // what the velocity fns expect; the distance is already in m.
                        let d_um = s.fiber_diameter_um;
                        let v_myel =
                            valenx_neuro::myelinated_conduction_velocity(d_um, s.k_myel_m_per_s_per_um);
                        let v_unmyel = valenx_neuro::unmyelinated_conduction_velocity(
                            d_um,
                            s.k_unmyel_m_per_s_per_sqrt_um,
                        );
                        let d_star = valenx_neuro::myelination_crossover_diameter(
                            s.k_myel_m_per_s_per_um,
                            s.k_unmyel_m_per_s_per_sqrt_um,
                        );
                        let delay_myel_ms =
                            valenx_neuro::conduction_delay_s(s.conduction_distance_m, v_myel) * 1000.0;
                        let delay_unmyel_ms =
                            valenx_neuro::conduction_delay_s(s.conduction_distance_m, v_unmyel)
                                * 1000.0;
                        ui.label(
                            egui::RichText::new(format!(
                                "myelinated {v_myel:.1} m/s · delay {delay_myel_ms:.2} ms\nunmyelinated {v_unmyel:.2} m/s · delay {delay_unmyel_ms:.2} ms\ncrossover diameter d* {d_star:.2} µm"
                            ))
                            .monospace()
                            .small(),
                        );
                    }

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
                        if !r.eap_uv.is_empty() {
                            let lo = r.eap_uv.iter().cloned().fold(f64::INFINITY, f64::min);
                            let hi = r.eap_uv.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                            ui.label(
                                egui::RichText::new(format!(
                                    "Recorded spike (EAP): {lo:.0} … {hi:.0} µV (biphasic)"
                                ))
                                .strong(),
                            );
                            Plot::new("neuro_eap").height(130.0).show(ui, |pui| {
                                let pts: Vec<[f64; 2]> = r
                                    .eap_uv
                                    .iter()
                                    .enumerate()
                                    .map(|(i, &v)| [i as f64 * r.eap_dt_ms, v])
                                    .collect();
                                pui.line(Line::new(PlotPoints::from(pts)).name("φe (µV)"));
                            });
                        }
                        ui.label(
                            egui::RichText::new("Cross-section (● electrode, — axons)")
                                .weak()
                                .small(),
                        );
                        draw_schematic(ui, r);
                    }
                });
        }, );
    if close {
        app.show_neuro_workbench = false;
    }
}

/// A 2-D cross-section: the electrode at top, axons at their depths, coloured
/// green when recruited.
fn draw_schematic(ui: &mut egui::Ui, r: &NeuroResults) {
    let (resp, painter) = ui.allocate_painter(
        egui::vec2(ui.available_width(), 120.0),
        egui::Sense::hover(),
    );
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
        assert!(!r.eap_uv.is_empty(), "EAP recorded");
        let lo = r.eap_uv.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = r.eap_uv.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            lo < 0.0 && hi > 0.0,
            "recorded EAP must be biphasic: {lo}..{hi}"
        );
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

    #[test]
    fn charge_safety_tracks_the_shannon_limit() {
        let mut s = NeuroWorkbenchState::default();
        let base = charge_safety(&s);
        assert!(base.q_uc > 0.0 && base.area_cm2 > 0.0 && base.k_value.is_finite());
        // Charge per phase scales linearly with current; half the max-safe
        // charge is safe, double is unsafe.
        let pw_factor = s.pulse_width_us * 1.0e-6; // µC per µA at this pulse width
        s.electrode_ua = 0.5 * base.max_safe_q_uc / pw_factor;
        assert!(
            charge_safety(&s).safe,
            "half the max-safe charge must be safe"
        );
        s.electrode_ua = 2.0 * base.max_safe_q_uc / pw_factor;
        assert!(
            !charge_safety(&s).safe,
            "double the max-safe charge must be unsafe"
        );
    }

    #[test]
    fn current_distance_readout_tracks_geometry() {
        let mut s = NeuroWorkbenchState {
            cd_i0_ua: 2.0,
            cd_k_ua_per_um2: 0.002,
            depth_mm: 0.5, // 500 µm
            ..Default::default()
        };
        // I_th(500 µm) = 2 + 0.002·500² = 502 µA.
        assert!((current_distance_readout(&s).threshold_at_depth_ua - 502.0).abs() < 1e-6);
        // 600 µA reaches the nearest axon; 400 µA does not.
        s.electrode_ua = -600.0;
        assert!(current_distance_readout(&s).reaches_nearest);
        s.electrode_ua = -400.0;
        assert!(!current_distance_readout(&s).reaches_nearest);
        // At exactly the threshold current, the activation radius == the depth.
        s.electrode_ua = -502.0;
        assert!((current_distance_readout(&s).activation_radius_um - 500.0).abs() < 1e-6);
    }
}
