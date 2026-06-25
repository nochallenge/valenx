//! The right-side **Sensors Workbench** panel — a native front-end over the
//! in-house `valenx-sensors` crate (LiDAR + radar ray-cast and measurement).
//!
//! Mirrors the other workbenches (`rotor_workbench`, `fem_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_sensors_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"sensors"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! The user picks a **sensor type** (LiDAR or Radar), edits a small scene
//! (a ground plane plus a sphere with adjustable position and radius), tunes
//! the sensor parameters, then clicks **Scan / Measure**. Results render in a
//! 2-D egui-painter visualisation:
//!
//! * **LiDAR** — top-down ray plot (hit points vs misses) + range stats.
//! * **Radar** — detection readout (range, SNR, Doppler, detected flag) + a
//!   simple range-vs-SNR bar.
//!
//! Honesty: the model is research / educational-grade. The LiDAR uses a
//! uniform beam grid; the radar is a simple monostatic link-budget
//! point-target model with additive noise. No sensor-fusion or tracking.
//! Every error from `valenx-sensors` surfaces verbatim — the workbench never
//! invents a number.

use eframe::egui;
use valenx_sensors::{
    AntennaPattern, Beam, Lidar, LidarConfig, Plane, Radar, RadarConfig, RadarReturn, Rcs, Scene,
    Sphere,
};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Sensor selector
// ---------------------------------------------------------------------------

/// Which sensor type the workbench is currently configured for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum SensorKind {
    #[default]
    Lidar,
    Radar,
}

// ---------------------------------------------------------------------------
// LiDAR parameters
// ---------------------------------------------------------------------------

/// Editable LiDAR configuration shown in the workbench.
#[derive(Clone, Debug)]
pub struct LidarParams {
    /// Horizontal field of view (degrees).
    pub h_fov_deg: f64,
    /// Vertical field of view (degrees).
    pub v_fov_deg: f64,
    /// Number of azimuth beams.
    pub azimuth_steps: usize,
    /// Number of elevation beams.
    pub elevation_steps: usize,
    /// Maximum detection range (m).
    pub max_range_m: f64,
}

impl Default for LidarParams {
    fn default() -> Self {
        Self {
            h_fov_deg: 90.0,
            v_fov_deg: 30.0,
            azimuth_steps: 32,
            elevation_steps: 8,
            max_range_m: 50.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Radar parameters
// ---------------------------------------------------------------------------

/// Editable radar configuration shown in the workbench.
#[derive(Clone, Debug)]
pub struct RadarParams {
    /// Transmit power (W).
    pub tx_power_w: f64,
    /// Carrier frequency (GHz) — converted to Hz for the API.
    pub freq_ghz: f64,
    /// Target position X (m) (Y=0, Z=0 scene axis).
    pub target_x_m: f64,
    /// Target radial velocity (m/s).
    pub target_vr_m_s: f64,
    /// RCS sphere radius (m) — uses the `Rcs::Sphere` model.
    pub rcs_sphere_radius_m: f64,
}

impl Default for RadarParams {
    fn default() -> Self {
        Self {
            tx_power_w: 1000.0,
            freq_ghz: 10.0,
            target_x_m: 25.0,
            target_vr_m_s: -10.0,
            rcs_sphere_radius_m: 0.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Scene parameters (shared between sensor types)
// ---------------------------------------------------------------------------

/// The small editable scene: a ground plane at `plane_y` plus one sphere.
#[derive(Clone, Debug)]
pub struct SceneParams {
    /// Sphere centre X (m).
    pub sphere_x_m: f64,
    /// Sphere centre Z (m, depth axis for the top-down LiDAR plot).
    pub sphere_z_m: f64,
    /// Sphere radius (m).
    pub sphere_radius_m: f64,
}

impl Default for SceneParams {
    fn default() -> Self {
        Self {
            sphere_x_m: 0.0,
            sphere_z_m: 15.0,
            sphere_radius_m: 2.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Last computed LiDAR scan (beam geometry cached for the painter).
#[derive(Default, Clone)]
pub struct LidarResult {
    pub scan: Vec<Beam>,
    pub num_returns: usize,
    pub min_range_m: f64,
    pub max_range_m: f64,
    pub mean_range_m: f64,
}

/// Last computed radar measurement.
#[derive(Default, Clone)]
pub struct RadarResult {
    /// `None` when sensor did not detect a target.
    pub ret: Option<RadarReturn>,
    pub max_range_m: f64,
}

/// Persistent state for the Sensors workbench.
pub struct SensorsWorkbenchState {
    /// Active sensor type.
    pub sensor_kind: SensorKind,
    /// LiDAR parameters.
    pub lidar: LidarParams,
    /// Radar parameters.
    pub radar: RadarParams,
    /// Scene parameters (sphere + ground plane).
    pub scene: SceneParams,
    /// Last LiDAR result (populated after a successful scan).
    pub lidar_result: Option<LidarResult>,
    /// Last radar result (populated after a successful measurement).
    pub radar_result: Option<RadarResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl Default for SensorsWorkbenchState {
    fn default() -> Self {
        Self {
            sensor_kind: SensorKind::Lidar,
            lidar: LidarParams::default(),
            radar: RadarParams::default(),
            scene: SceneParams::default(),
            lidar_result: None,
            radar_result: None,
            status: String::new(),
        }
    }
}

impl SensorsWorkbenchState {
    /// Build the shared [`Scene`] from the current scene parameters.
    ///
    /// Returns `Err` (instead of panicking) when a geometry constructor rejects
    /// the current parameter values — most commonly when the user has dragged
    /// `sphere_radius_m` to ≤ 0 or to NaN.
    fn build_scene(&self) -> Result<Scene, String> {
        use nalgebra::Vector3;
        let mut scene = Scene::new();
        // Ground plane at Y = -1.0, normal pointing up.
        scene.push_plane(
            Plane::new(Vector3::new(0.0, -1.0, 0.0), Vector3::new(0.0, 1.0, 0.0))
                .map_err(|e| e.to_string())?,
        );
        // Sphere at the user's position.
        scene.push_sphere(
            Sphere::new(
                Vector3::new(self.scene.sphere_x_m, 0.0, self.scene.sphere_z_m),
                self.scene.sphere_radius_m,
            )
            .map_err(|e| e.to_string())?,
        );
        Ok(scene)
    }

    /// Run a LiDAR scan. Returns the populated [`LidarResult`] or an error string.
    pub fn scan_lidar(&self) -> Result<LidarResult, String> {
        use nalgebra::{UnitQuaternion, Vector3};
        let cfg = LidarConfig {
            azimuth_steps: self.lidar.azimuth_steps,
            elevation_steps: self.lidar.elevation_steps,
            h_fov: self.lidar.h_fov_deg.to_radians(),
            v_fov: self.lidar.v_fov_deg.to_radians(),
            min_range: 0.2,
            max_range: self.lidar.max_range_m,
            range_noise_std: 0.01,
        };
        let mut sensor = Lidar::new(cfg, 42).map_err(|e| e.to_string())?;
        let scene = self.build_scene()?;
        let scan = sensor.scan(&scene, Vector3::zeros(), UnitQuaternion::identity());

        let ranges: Vec<f64> = scan.beams.iter().filter_map(|b| b.range).collect();
        let num_returns = ranges.len();
        let (min_r, max_r, mean_r) = if ranges.is_empty() {
            (0.0, 0.0, 0.0)
        } else {
            let mn = ranges.iter().cloned().fold(f64::INFINITY, f64::min);
            let mx = ranges.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let avg = ranges.iter().sum::<f64>() / ranges.len() as f64;
            (mn, mx, avg)
        };

        Ok(LidarResult {
            scan: scan.beams,
            num_returns,
            min_range_m: min_r,
            max_range_m: max_r,
            mean_range_m: mean_r,
        })
    }

    /// Run a radar measurement. Returns the populated [`RadarResult`] or an error string.
    pub fn measure_radar(&self) -> Result<RadarResult, String> {
        use nalgebra::Vector3;
        let wavelength = 3e8 / (self.radar.freq_ghz * 1e9);
        let cfg = RadarConfig {
            tx_power_w: self.radar.tx_power_w,
            antenna: AntennaPattern::Gaussian {
                peak_gain: 1000.0, // 30 dBi
                beamwidth: 0.05,   // 3°
            },
            wavelength,
            bandwidth_hz: 1e6,
            noise_figure_db: 3.0,
            detection_threshold_db: 13.0,
            noise_std_db: 0.5,
        };
        let mut sensor = Radar::new(cfg, 99).map_err(|e| e.to_string())?;
        let rcs_sigma = Rcs::Sphere {
            radius: self.radar.rcs_sphere_radius_m,
        }
        .sigma()
        .map_err(|e| e.to_string())?;
        let max_range = sensor
            .max_detection_range(rcs_sigma)
            .map_err(|e| e.to_string())?;
        let target_pos = Vector3::new(self.radar.target_x_m, 0.0, 0.0);
        let target_vel = Vector3::new(self.radar.target_vr_m_s, 0.0, 0.0);
        let ret = sensor
            .measure(target_pos, target_vel, rcs_sigma)
            .map_err(|e| e.to_string())?;
        Ok(RadarResult {
            ret,
            max_range_m: max_range,
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`.
    /// Covers the always-visible scene controls plus BOTH sensor families
    /// (LiDAR + radar) so an agent can set either regardless of which sensor
    /// kind is currently selected; each string matches the form caption exactly.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            // -- Scene --
            "sphere X (m)",
            "sphere Z (m)",
            "sphere radius (m)",
            // -- LiDAR --
            "H field-of-view (°)",
            "V field-of-view (°)",
            "azimuth beams",
            "elevation beams",
            "max range (m)",
            // -- Radar --
            "Tx power (W)",
            "frequency (GHz)",
            "target X (m)",
            "target radial vel (m/s)",
            "RCS sphere radius (m)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` — never a panic. Ranges
    /// mirror the form's `DragValue` clamps exactly. The caption routes the
    /// value to the right sub-struct (`scene` / `lidar` / `radar`).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let ranged = |v: f64, lo: f64, hi: f64, what: &str| -> Result<f64, String> {
            if v.is_finite() && (lo..=hi).contains(&v) {
                Ok(v)
            } else {
                Err(format!("{what} must be in {lo}..={hi}, got {v}"))
            }
        };
        let ranged_int = |value: &crate::agent_commands::AgentValue,
                          lo: i64,
                          hi: i64,
                          what: &str|
         -> Result<usize, String> {
            let n = value.as_i64()?;
            if (lo..=hi).contains(&n) {
                Ok(n as usize)
            } else {
                Err(format!("{what} must be in {lo}..={hi}, got {n}"))
            }
        };
        match name {
            // -- Scene --
            "sphere X (m)" => {
                self.scene.sphere_x_m = ranged(value.as_f64()?, -100.0, 100.0, "sphere X")?
            }
            "sphere Z (m)" => {
                self.scene.sphere_z_m = ranged(value.as_f64()?, 0.1, 200.0, "sphere Z")?
            }
            "sphere radius (m)" => {
                self.scene.sphere_radius_m = ranged(value.as_f64()?, 0.01, 50.0, "sphere radius")?
            }
            // -- LiDAR --
            "H field-of-view (°)" => {
                self.lidar.h_fov_deg = ranged(value.as_f64()?, 1.0, 360.0, "H field-of-view")?
            }
            "V field-of-view (°)" => {
                self.lidar.v_fov_deg = ranged(value.as_f64()?, 1.0, 90.0, "V field-of-view")?
            }
            "azimuth beams" => {
                self.lidar.azimuth_steps = ranged_int(value, 4, 512, "azimuth beams")?
            }
            "elevation beams" => {
                self.lidar.elevation_steps = ranged_int(value, 1, 64, "elevation beams")?
            }
            "max range (m)" => {
                self.lidar.max_range_m = ranged(value.as_f64()?, 1.0, 500.0, "max range")?
            }
            // -- Radar --
            "Tx power (W)" => {
                self.radar.tx_power_w = ranged(value.as_f64()?, 1.0, 1e6, "Tx power")?
            }
            "frequency (GHz)" => {
                self.radar.freq_ghz = ranged(value.as_f64()?, 0.1, 100.0, "frequency")?
            }
            "target X (m)" => {
                self.radar.target_x_m = ranged(value.as_f64()?, 0.1, 500.0, "target X")?
            }
            "target radial vel (m/s)" => {
                self.radar.target_vr_m_s =
                    ranged(value.as_f64()?, -300.0, 300.0, "target radial vel")?
            }
            "RCS sphere radius (m)" => {
                self.radar.rcs_sphere_radius_m =
                    ranged(value.as_f64()?, 0.01, 20.0, "RCS sphere radius")?
            }
            other => return Err(format!("unknown Sensors control: {other:?}")),
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Sensors workbench (a no-op unless toggled on via View → Sensors).
/// Mirrors [`crate::rotor_workbench::draw_rotor_workbench`].
pub fn draw_sensors_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_sensors_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_sensors_workbench",
        "Sensors (LiDAR / Radar)",
        sensors_workbench_body,
    );
    if close {
        app.show_sensors_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn sensors_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new("LiDAR + radar sensor suite · valenx-sensors")
            .weak()
            .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.sensors;

        // --- Sensor selector -------------------------------------------------
        ui.label(egui::RichText::new("Sensor type").strong());
        ui.horizontal(|ui| {
            ui.radio_value(&mut s.sensor_kind, SensorKind::Lidar, "LiDAR");
            ui.radio_value(&mut s.sensor_kind, SensorKind::Radar, "Radar");
        });
        ui.add_space(6.0);

        // --- Scene params ----------------------------------------------------
        ui.label(egui::RichText::new("Scene").strong());
        egui::Grid::new("sensors_scene")
            .num_columns(2)
            .show(ui, |ui| {
                let lbl = ui.label("sphere X (m)");
                ui.add(
                    egui::DragValue::new(&mut s.scene.sphere_x_m)
                        .speed(0.5)
                        .range(-100.0..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Sphere centre X position (m)");
                ui.end_row();

                let lbl = ui.label("sphere Z (m)");
                ui.add(
                    egui::DragValue::new(&mut s.scene.sphere_z_m)
                        .speed(0.5)
                        .range(0.1..=200.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Sphere centre Z position / depth (m)");
                ui.end_row();

                let lbl = ui.label("sphere radius (m)");
                ui.add(
                    egui::DragValue::new(&mut s.scene.sphere_radius_m)
                        .speed(0.1)
                        .range(0.01..=50.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Sphere radius (m)");
                ui.end_row();
            });

        ui.add_space(6.0);

        // --- Sensor-specific params ------------------------------------------
        match s.sensor_kind {
            SensorKind::Lidar => {
                ui.label(egui::RichText::new("LiDAR parameters").strong());
                egui::Grid::new("sensors_lidar_params")
                    .num_columns(2)
                    .show(ui, |ui| {
                        let lbl = ui.label("H field-of-view (°)");
                        ui.add(
                            egui::DragValue::new(&mut s.lidar.h_fov_deg)
                                .speed(1.0)
                                .range(1.0..=360.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Horizontal field of view (degrees)");
                        ui.end_row();

                        let lbl = ui.label("V field-of-view (°)");
                        ui.add(
                            egui::DragValue::new(&mut s.lidar.v_fov_deg)
                                .speed(1.0)
                                .range(1.0..=90.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Vertical field of view (degrees)");
                        ui.end_row();

                        let lbl = ui.label("azimuth beams");
                        ui.add(
                            egui::DragValue::new(&mut s.lidar.azimuth_steps)
                                .speed(1)
                                .range(4..=512),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Number of azimuth (horizontal) beam steps");
                        ui.end_row();

                        let lbl = ui.label("elevation beams");
                        ui.add(
                            egui::DragValue::new(&mut s.lidar.elevation_steps)
                                .speed(1)
                                .range(1..=64),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Number of elevation (vertical) beam steps");
                        ui.end_row();

                        let lbl = ui.label("max range (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.lidar.max_range_m)
                                .speed(1.0)
                                .range(1.0..=500.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Maximum LiDAR detection range (m)");
                        ui.end_row();
                    });
            }
            SensorKind::Radar => {
                ui.label(egui::RichText::new("Radar parameters").strong());
                egui::Grid::new("sensors_radar_params")
                    .num_columns(2)
                    .show(ui, |ui| {
                        let lbl = ui.label("Tx power (W)");
                        ui.add(
                            egui::DragValue::new(&mut s.radar.tx_power_w)
                                .speed(10.0)
                                .range(1.0..=1e6),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Transmit power (watts)");
                        ui.end_row();

                        let lbl = ui.label("frequency (GHz)");
                        ui.add(
                            egui::DragValue::new(&mut s.radar.freq_ghz)
                                .speed(0.5)
                                .range(0.1..=100.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Carrier frequency (GHz)");
                        ui.end_row();

                        let lbl = ui.label("target X (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.radar.target_x_m)
                                .speed(1.0)
                                .range(0.1..=500.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Target range along X axis (m)");
                        ui.end_row();

                        let lbl = ui.label("target radial vel (m/s)");
                        ui.add(
                            egui::DragValue::new(&mut s.radar.target_vr_m_s)
                                .speed(0.5)
                                .range(-300.0..=300.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text("Target radial velocity (m/s); negative = approaching");
                        ui.end_row();

                        let lbl = ui.label("RCS sphere radius (m)");
                        ui.add(
                            egui::DragValue::new(&mut s.radar.rcs_sphere_radius_m)
                                .speed(0.05)
                                .range(0.01..=20.0),
                        )
                        .labelled_by(lbl.id)
                        .on_hover_text(
                            "Sphere radius used for the radar cross-section (Mie) model (m)",
                        );
                        ui.end_row();
                    });
            }
        }

        ui.add_space(6.0);
        let btn_label = match s.sensor_kind {
            SensorKind::Lidar => "Scan",
            SensorKind::Radar => "Measure",
        };
        if ui
            .button(egui::RichText::new(btn_label).strong())
            .on_hover_text("Run the sensor model against the current scene.")
            .clicked()
        {
            do_run = true;
        }
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        let s = &mut app.sensors;
        match s.sensor_kind {
            SensorKind::Lidar => match s.scan_lidar() {
                Ok(r) => {
                    s.status = format!(
                        "✔ {} hits / {} beams · range [{:.2}, {:.2}] m · mean {:.2} m",
                        r.num_returns,
                        r.scan.len(),
                        r.min_range_m,
                        r.max_range_m,
                        r.mean_range_m,
                    );
                    s.lidar_result = Some(r);
                    s.radar_result = None;
                }
                Err(e) => {
                    s.status = format!("⚠ {e}");
                    s.lidar_result = None;
                }
            },
            SensorKind::Radar => match s.measure_radar() {
                Ok(r) => {
                    if let Some(ret) = &r.ret {
                        s.status = format!(
                            "✔ detected · range {:.2} m · SNR {:.1} dB · Doppler {:.1} Hz",
                            ret.range, ret.snr_db, ret.doppler_hz,
                        );
                    } else {
                        s.status = format!("○ no detection (max range ≈ {:.1} m)", r.max_range_m);
                    }
                    s.radar_result = Some(r);
                    s.lidar_result = None;
                }
                Err(e) => {
                    s.status = format!("⚠ {e}");
                    s.radar_result = None;
                }
            },
        }
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.sensors;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('⚠') {
            egui::Color32::from_rgb(220, 120, 60)
        } else if s.status.starts_with('○') {
            egui::Color32::from_rgb(200, 160, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();

    match s.sensor_kind {
        SensorKind::Lidar => draw_lidar_viz(s, ui),
        SensorKind::Radar => draw_radar_viz(s, ui),
    }
}

// ---------------------------------------------------------------------------
// LiDAR visualisation (top-down 2-D painter)
// ---------------------------------------------------------------------------

fn draw_lidar_viz(s: &SensorsWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("LiDAR top-down view").strong());
    ui.label(
        egui::RichText::new("green = hit · grey = miss · axes: X (right) / Z (down)")
            .weak()
            .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(400.0), 280.0),
        egui::Sense::hover(),
    );

    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(20, 22, 28));

    let Some(res) = &s.lidar_result else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "press Scan to visualise",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(110),
        );
        return;
    };

    // Map scene coords (X horizontal, Z depth) onto the painter rect.
    // Scene range: X ∈ [-max_range, max_range], Z ∈ [0, max_range].
    let max_r = s.lidar.max_range_m as f32;
    let margin = 10.0_f32;
    let inner = rect.shrink(margin);

    let scene_to_px = |x: f32, z: f32| -> egui::Pos2 {
        let nx = (x + max_r) / (2.0 * max_r); // 0..1
        let nz = z / max_r;
        egui::Pos2::new(
            inner.left() + nx * inner.width(),
            inner.top() + nz * inner.height(),
        )
    };

    // Origin dot.
    painter.circle_filled(
        scene_to_px(0.0, 0.0),
        4.0,
        egui::Color32::from_rgb(220, 220, 80),
    );

    // Draw sphere outline for reference.
    let sc = scene_to_px(s.scene.sphere_x_m as f32, s.scene.sphere_z_m as f32);
    let sr = (s.scene.sphere_radius_m as f32 / (2.0 * max_r)) * inner.width();
    painter.circle_stroke(
        sc,
        sr.max(2.0),
        egui::Stroke::new(1.0, egui::Color32::from_rgb(100, 160, 240)),
    );

    // Draw beams: only the azimuth projection (take elevation_step 0 per az, or all beams).
    for beam in &res.scan {
        // direction projected onto X/Z plane
        let dx = beam.direction.x as f32;
        let dz = beam.direction.z as f32;
        if let Some(r) = beam.range {
            let r = r as f32;
            let hit = scene_to_px(dx * r, dz * r);
            painter.line_segment(
                [scene_to_px(0.0, 0.0), hit],
                egui::Stroke::new(0.5, egui::Color32::from_rgba_unmultiplied(60, 200, 80, 80)),
            );
            painter.circle_filled(hit, 2.0, egui::Color32::from_rgb(60, 220, 80));
        } else {
            let far = scene_to_px(dx * max_r * 0.95, dz * max_r * 0.95);
            painter.line_segment(
                [scene_to_px(0.0, 0.0), far],
                egui::Stroke::new(
                    0.3,
                    egui::Color32::from_rgba_unmultiplied(120, 120, 120, 40),
                ),
            );
        }
    }

    // Stats
    ui.add_space(4.0);
    egui::Grid::new("sensors_lidar_stats")
        .num_columns(2)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "hits", format!("{}", res.num_returns));
            row(ui, "beams", format!("{}", res.scan.len()));
            row(ui, "min range (m)", format!("{:.3}", res.min_range_m));
            row(ui, "max range (m)", format!("{:.3}", res.max_range_m));
            row(ui, "mean range (m)", format!("{:.3}", res.mean_range_m));
        });
}

// ---------------------------------------------------------------------------
// Radar visualisation (detection readout + range-SNR indicator)
// ---------------------------------------------------------------------------

fn draw_radar_viz(s: &SensorsWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Radar measurement").strong());

    let Some(res) = &s.radar_result else {
        ui.label(
            egui::RichText::new("press Measure to compute a radar return")
                .weak()
                .small(),
        );
        return;
    };

    // Detection readout.
    egui::Grid::new("sensors_radar_readout")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(
                ui,
                "detected",
                if res.ret.is_some() {
                    "YES".to_string()
                } else {
                    "NO".to_string()
                },
            );
            row(
                ui,
                "max detection range (m)",
                format!("{:.1}", res.max_range_m),
            );
            if let Some(ret) = &res.ret {
                row(ui, "measured range (m)", format!("{:.3}", ret.range));
                row(ui, "SNR (dB)", format!("{:.2}", ret.snr_db));
                row(ui, "Doppler (Hz)", format!("{:.2}", ret.doppler_hz));
            }
        });

    ui.add_space(8.0);

    // Simple range-vs-SNR indicator bar.
    ui.label(egui::RichText::new("SNR vs. max-range fraction").strong());
    let available_x = ui.available_width().min(400.0);
    let (bar_rect, _) = ui.allocate_exact_size(egui::vec2(available_x, 28.0), egui::Sense::hover());
    let painter = ui.painter_at(bar_rect);
    painter.rect_filled(bar_rect, 4.0, egui::Color32::from_rgb(30, 32, 40));

    if let Some(ret) = &res.ret {
        let fraction = (ret.range / res.max_range_m.max(1.0)).clamp(0.0, 1.0) as f32;
        let snr_norm = ((ret.snr_db + 10.0) / 50.0).clamp(0.0, 1.0) as f32;
        // Range bar (blue): shows how far the target is relative to max range.
        let range_w = bar_rect.width() * fraction;
        painter.rect_filled(
            egui::Rect::from_min_size(bar_rect.left_top(), egui::vec2(range_w, bar_rect.height())),
            2.0,
            egui::Color32::from_rgb(60, 120, 220),
        );
        // SNR bar (green, inner).
        let snr_w = bar_rect.width() * snr_norm;
        painter.rect_filled(
            egui::Rect::from_min_size(
                bar_rect.left_top() + egui::vec2(0.0, bar_rect.height() * 0.6),
                egui::vec2(snr_w, bar_rect.height() * 0.4),
            ),
            0.0,
            egui::Color32::from_rgb(60, 200, 90),
        );
    }
    painter.text(
        bar_rect.left_center() + egui::vec2(6.0, 0.0),
        egui::Align2::LEFT_CENTER,
        "range",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(200),
    );
    painter.text(
        bar_rect.right_center() + egui::vec2(-6.0, 0.0),
        egui::Align2::RIGHT_CENTER,
        "max range",
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(160),
    );
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = SensorsWorkbenchState::default();
        // Scene float + LiDAR integer both land in the right sub-struct.
        s.agent_set("sphere radius (m)", &AgentValue::Float(2.0))
            .unwrap();
        assert_eq!(s.scene.sphere_radius_m, 2.0);
        s.agent_set("azimuth beams", &AgentValue::Int(128)).unwrap();
        assert_eq!(s.lidar.azimuth_steps, 128);
        // A radar caption routes to the radar sub-struct.
        s.agent_set("Tx power (W)", &AgentValue::Float(500.0))
            .unwrap();
        assert_eq!(s.radar.tx_power_w, 500.0);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into a numeric field) -> Err.
        assert!(s
            .agent_set("sphere radius (m)", &AgentValue::Str("big".into()))
            .is_err());
        // Out-of-range (azimuth beams > 512) -> Err, field untouched.
        assert!(s
            .agent_set("azimuth beams", &AgentValue::Int(9999))
            .is_err());
        assert_eq!(
            s.lidar.azimuth_steps, 128,
            "rejected set leaves field untouched"
        );
    }

    #[test]
    fn default_lidar_scan_returns_finite_ranges() {
        let s = SensorsWorkbenchState::default();
        let res = s.scan_lidar().expect("default LiDAR scan should succeed");
        assert!(
            res.num_returns > 0,
            "default scene has a sphere in range — expect hits"
        );
        assert!(res.min_range_m.is_finite() && res.min_range_m > 0.0);
        assert!(res.max_range_m >= res.min_range_m);
        assert!(res.mean_range_m >= res.min_range_m && res.mean_range_m <= res.max_range_m);
    }

    #[test]
    fn lidar_scan_beam_count_matches_config() {
        let s = SensorsWorkbenchState::default();
        let res = s.scan_lidar().expect("scan should succeed");
        let expected = s.lidar.azimuth_steps * s.lidar.elevation_steps;
        assert_eq!(res.scan.len(), expected, "beam count must match config");
    }

    #[test]
    fn default_radar_detects_close_target() {
        let s = SensorsWorkbenchState::default();
        // Default target at 25 m should be well within range.
        let res = s.measure_radar().expect("radar measure should succeed");
        assert!(
            res.ret.is_some(),
            "default target at 25 m should be detected (max range ≈ {:.1} m)",
            res.max_range_m
        );
    }

    #[test]
    fn radar_max_range_is_finite_and_positive() {
        // max_detection_range must be finite and positive for the default config.
        let s = SensorsWorkbenchState::default();
        let res = s.measure_radar().expect("measure should succeed");
        assert!(
            res.max_range_m.is_finite() && res.max_range_m > 0.0,
            "max detection range must be finite and positive, got {}",
            res.max_range_m
        );
    }

    #[test]
    fn bad_radar_config_zero_power_errors() {
        let mut s = SensorsWorkbenchState::default();
        s.radar.tx_power_w = 0.0; // invalid: must be > 0
        let res = s.measure_radar();
        assert!(res.is_err(), "zero tx_power_w should produce an error");
    }

    #[test]
    fn bad_lidar_config_zero_azimuth_steps_errors() {
        let mut s = SensorsWorkbenchState::default();
        s.lidar.azimuth_steps = 0; // invalid
        let res = s.scan_lidar();
        // valenx-sensors should reject zero-step config
        assert!(res.is_err(), "zero azimuth steps should produce an error");
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
            draw_sensors_workbench(app, ctx);
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
        assert!(!app.show_sensors_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_sensors_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_lidar_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        // Pre-populate a result so the painter branch executes.
        let res = app.sensors.scan_lidar().expect("scan should succeed");
        app.sensors.lidar_result = Some(res);
        app.sensors.status = "✔ test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_radar_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        app.sensors.sensor_kind = SensorKind::Radar;
        let res = app.sensors.measure_radar().expect("measure should succeed");
        app.sensors.radar_result = Some(res);
        app.sensors.status = "✔ radar test".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        app.sensors.status = "⚠ simulated error for testing".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();

        // 3 scene params + 5 LiDAR params (default sensor kind) = 8 minimum.
        assert!(
            spin_buttons.len() >= 8,
            "expected at least 8 numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        for caption in [
            "sphere X (m)",
            "sphere Z (m)",
            "sphere radius (m)",
            "max range (m)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Scan button must be named and invokable.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Scan"))
            }),
            "the Scan button must be a named, invokable node"
        );
    }

    #[test]
    fn radar_numeric_controls_are_named() {
        let mut app = ValenxApp::default();
        app.show_sensors_workbench = true;
        app.sensors.sensor_kind = SensorKind::Radar;
        let nodes = draw_and_collect_nodes(&mut app);

        for caption in ["Tx power (W)", "frequency (GHz)", "target X (m)"] {
            assert!(
                has_named_node(&nodes, caption),
                "radar caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // Measure button must be present.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Measure"))
            }),
            "the Measure button must be a named, invokable node"
        );
    }

    #[test]
    fn invalid_sphere_radius_returns_err_not_panic() {
        // When the user drags sphere_radius_m to <= 0, build_scene must return
        // Err and scan_lidar must propagate it — NOT panic.
        let mut state = SensorsWorkbenchState::default();
        state.scene.sphere_radius_m = 0.0;
        assert!(
            state.scan_lidar().is_err(),
            "sphere_radius_m = 0 must produce Err, not panic"
        );
        state.scene.sphere_radius_m = -1.5;
        assert!(
            state.scan_lidar().is_err(),
            "sphere_radius_m < 0 must produce Err, not panic"
        );
    }

    #[test]
    fn agent_bridge_sensors_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge (agent_commands::apply)
        // uses for `OpenWorkbench { id: "sensors" }`:
        //   1. TabKind::from_id("sensors") → Some(TabKind::Sensors)
        //   2. set_workbench_flag(app, "sensors", true) → show_sensors_workbench = true
        // The agent_commands::apply fn is module-private, so we test its
        // constituent public pieces; open_workbench_switches_the_active_tab_kind
        // in agent_commands.rs already covers the composed path for other kinds.
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("sensors"),
            Some(TabKind::Sensors),
            "\"sensors\" must resolve to TabKind::Sensors"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("SENSORS"), Some(TabKind::Sensors));
        assert_eq!(TabKind::from_id("  sensors  "), Some(TabKind::Sensors));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_sensors_workbench);
        set_workbench_flag(&mut app, "sensors", true);
        assert!(
            app.show_sensors_workbench,
            "set_workbench_flag(\"sensors\", true) must set show_sensors_workbench"
        );
        set_workbench_flag(&mut app, "sensors", false);
        assert!(!app.show_sensors_workbench);
    }
}
