//! The right-side **Car** workbench — a design → simulate panel over
//! [`valenx_vehicle`]. Pick a preset or tune the performance-relevant
//! parameters (mass, power, drivetrain, grip, aero) and the panel
//! reactively reports top speed, the 0-100 / 0-200 km/h sprints, the
//! 100-0 braking distance, and steady-state skidpad grip, plus a
//! full-throttle speed-vs-time launch curve.
//!
//! Mirrors the other single-file workbenches (`rocket_workbench`,
//! `engine_workbench`): a resizable [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_car_workbench`, toggled from the View menu or
//! opened as a project tab.
//!
//! Honest scope: research / preliminary-design grade — a point-mass
//! longitudinal + friction-circle model (power- or traction-limited
//! tractive force with longitudinal weight transfer and speed-dependent
//! downforce). It is not a full vehicle-dynamics or lap-time simulation.

use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::ValenxApp;
use valenx_vehicle::{hypercar, sports_car, Car, Drivetrain};

/// Metres-per-second to kilometres-per-hour.
const MS_TO_KMH: f64 = 3.6;
/// Skidpad corner radius the readout reports against (m).
const SKIDPAD_RADIUS_M: f64 = 30.0;

/// Computed performance for one [`Car`].
struct CarPerformance {
    top_speed_kmh: f64,
    sprint_0_100_s: f64,
    sprint_0_200_s: f64,
    braking_100_0_m: f64,
    skidpad_kmh: f64,
    skidpad_lat_g: f64,
    /// `[time_s, speed_kmh]` full-throttle standing-start launch curve.
    launch_curve: Vec<[f64; 2]>,
}

/// Run the full performance suite for `car`.
fn simulate(car: &Car) -> CarPerformance {
    let top = car.top_speed();
    let skid_v = car.skidpad_speed(SKIDPAD_RADIUS_M);
    CarPerformance {
        top_speed_kmh: top * MS_TO_KMH,
        sprint_0_100_s: car.accelerate_to(100.0 / MS_TO_KMH).time,
        sprint_0_200_s: car.accelerate_to(200.0 / MS_TO_KMH).time,
        braking_100_0_m: car.braking_distance(100.0 / MS_TO_KMH),
        skidpad_kmh: skid_v * MS_TO_KMH,
        skidpad_lat_g: car.max_lateral_g(skid_v),
        launch_curve: launch_curve(car, top),
    }
}

/// `[time_s, speed_kmh]` samples of a standing-start full-throttle run up
/// to (just below) top speed. Coarse `dt` — this feeds a plot, not a
/// timing figure (the sprint times use the crate's fine integrator).
fn launch_curve(car: &Car, top: f64) -> Vec<[f64; 2]> {
    let dt = 0.01;
    let cap = top * 0.999;
    let (mut v, mut t) = (0.0_f64, 0.0_f64);
    let mut curve = vec![[0.0, 0.0]];
    while v < cap && t < 60.0 {
        let a = car.acceleration_at(v);
        if a <= 1e-4 {
            break;
        }
        v += a * dt;
        t += dt;
        curve.push([t, v * MS_TO_KMH]);
    }
    curve
}

/// Persistent form + result state for the Car workbench.
pub struct CarWorkbenchState {
    /// The car being designed.
    car: Car,
    /// The car the cached [`Self::performance`] was computed for; drives
    /// the reactive recompute (recompute only when it differs from `car`).
    last: Option<Car>,
    /// Cached performance for `last`.
    performance: Option<CarPerformance>,
}

impl Default for CarWorkbenchState {
    fn default() -> Self {
        Self {
            car: sports_car(),
            last: None,
            performance: None,
        }
    }
}

/// Draw the Car workbench panel (right dock), gated on the show flag.
pub fn draw_car_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_car_workbench {
        return;
    }
    egui::SidePanel::right("valenx_car_workbench")
        .resizable(true)
        .default_width(380.0)
        .width_range(330.0..=620.0)
        .show(ctx, |ui| {
            if crate::workbench_ui::header(
                ui,
                "Car — design → simulate",
                "point-mass longitudinal + friction-circle · valenx-vehicle",
            ) {
                app.show_car_workbench = false;
            }
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    let s = &mut app.car;

                    ui.horizontal(|ui| {
                        ui.label("preset:");
                        if ui.button("Sports car").clicked() {
                            s.car = sports_car();
                        }
                        if ui.button("Hypercar").clicked() {
                            s.car = hypercar();
                        }
                    });
                    ui.separator();

                    ui.add(egui::Slider::new(&mut s.car.mass, 600.0..=3000.0).text("mass (kg)"));
                    let mut power_kw = s.car.peak_power / 1000.0;
                    if ui
                        .add(
                            egui::Slider::new(&mut power_kw, 50.0..=1500.0).text("peak power (kW)"),
                        )
                        .changed()
                    {
                        s.car.peak_power = power_kw * 1000.0;
                    }
                    ui.horizontal(|ui| {
                        ui.label("drivetrain:");
                        ui.radio_value(&mut s.car.drivetrain, Drivetrain::FrontWheel, "FWD");
                        ui.radio_value(&mut s.car.drivetrain, Drivetrain::RearWheel, "RWD");
                        ui.radio_value(&mut s.car.drivetrain, Drivetrain::AllWheel, "AWD");
                    });
                    ui.add(
                        egui::Slider::new(&mut s.car.tire_friction, 0.6..=2.0).text("tire grip μ"),
                    );
                    ui.add(
                        egui::Slider::new(&mut s.car.drag_coefficient, 0.20..=0.60).text("drag Cd"),
                    );
                    ui.add(
                        egui::Slider::new(&mut s.car.frontal_area, 1.5..=3.0)
                            .text("frontal area (m²)"),
                    );
                    ui.add(
                        egui::Slider::new(&mut s.car.downforce_cla, 0.0..=6.0)
                            .text("downforce Cl·A (m²)"),
                    );

                    // Reactive recompute, only when the design actually changed.
                    if s.last != Some(s.car) {
                        s.performance = Some(simulate(&s.car));
                        s.last = Some(s.car);
                    }

                    ui.separator();
                    if let Some(p) = &s.performance {
                        ui.label(
                            egui::RichText::new(format!(
                                "top speed    {:>6.1} km/h\n\
                                 0-100 km/h   {:>6.2} s\n\
                                 0-200 km/h   {:>6.2} s\n\
                                 100-0 brake  {:>6.1} m\n\
                                 skidpad 30 m {:>6.1} km/h   ({:.2} g)",
                                p.top_speed_kmh,
                                p.sprint_0_100_s,
                                p.sprint_0_200_s,
                                p.braking_100_0_m,
                                p.skidpad_kmh,
                                p.skidpad_lat_g,
                            ))
                            .monospace()
                            .small(),
                        );
                        if p.launch_curve.len() > 1 {
                            ui.label(
                                egui::RichText::new(
                                    "speed (km/h) vs time (s) — full-throttle launch",
                                )
                                .weak()
                                .small(),
                            );
                            Plot::new("car_launch_plot").height(200.0).show(ui, |pui| {
                                pui.line(
                                    Line::new(PlotPoints::from(p.launch_curve.clone()))
                                        .name("speed (km/h)"),
                                );
                            });
                        }
                    }
                });
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sports_car_performance_is_physically_sane() {
        let p = simulate(&sports_car());
        // ~400 hp RWD sports car: high top speed, brisk 0-100, real braking.
        assert!(
            p.top_speed_kmh > 200.0 && p.top_speed_kmh < 400.0,
            "top speed {} km/h",
            p.top_speed_kmh
        );
        assert!(
            p.sprint_0_100_s > 2.0 && p.sprint_0_100_s < 9.0,
            "0-100 {} s",
            p.sprint_0_100_s
        );
        assert!(
            p.sprint_0_200_s > p.sprint_0_100_s,
            "0-200 must be slower than 0-100"
        );
        assert!(
            p.braking_100_0_m > 20.0 && p.braking_100_0_m < 90.0,
            "braking {} m",
            p.braking_100_0_m
        );
        assert!(p.skidpad_lat_g > 0.8, "skidpad {} g", p.skidpad_lat_g);
        assert!(p.launch_curve.len() > 5, "launch curve has samples");
    }

    #[test]
    fn hypercar_out_accelerates_the_sports_car() {
        let sc = simulate(&sports_car());
        let hc = simulate(&hypercar());
        assert!(
            hc.sprint_0_100_s < sc.sprint_0_100_s,
            "hypercar 0-100 ({}) should beat sports car ({})",
            hc.sprint_0_100_s,
            sc.sprint_0_100_s
        );
        assert!(
            hc.top_speed_kmh > sc.top_speed_kmh,
            "hypercar should have a higher top speed"
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    fn draw(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_car_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_car_workbench);
        draw(&mut app);
        assert!(app.car.performance.is_none());
    }

    #[test]
    fn workbench_computes_and_draws_on_first_open() {
        let mut app = ValenxApp::default();
        app.show_car_workbench = true;
        draw(&mut app);
        assert!(
            app.car.performance.is_some(),
            "first draw computes performance"
        );
    }
}
