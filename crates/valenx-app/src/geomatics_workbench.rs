//! The right-side **Geomatics Workbench** panel — native geodesic
//! calculations over `valenx-geomatics`.
//!
//! Mirrors the springs / gears / CFD workbenches: a resizable
//! [`egui::SidePanel`] gated on
//! `crate::ValenxApp::show_geomatics_workbench`, toggled from the View
//! menu. The form takes three WGS84 points — A, B, and a query point P —
//! and the "Compute" button reports the great-circle (haversine) and
//! rhumb-line distances A→B, the initial and final bearings, and the
//! cross-track / along-track offsets of P from the path A→B, as a
//! monospace readout.

use eframe::egui;

use valenx_geomatics::{
    along_track_distance, cross_track_distance, final_bearing, haversine_distance, initial_bearing,
    rhumb_distance, LatLon,
};

use crate::ValenxApp;

/// Persistent form + result state for the Geomatics Workbench.
pub struct GeomaticsWorkbenchState {
    /// Point A latitude (decimal degrees).
    lat_a: f64,
    /// Point A longitude (decimal degrees).
    lon_a: f64,
    /// Point B latitude (decimal degrees).
    lat_b: f64,
    /// Point B longitude (decimal degrees).
    lon_b: f64,
    /// Query point P latitude (decimal degrees).
    lat_p: f64,
    /// Query point P longitude (decimal degrees).
    lon_p: f64,
    /// Formatted geodesic readout (empty until the first compute).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
}

impl Default for GeomaticsWorkbenchState {
    fn default() -> Self {
        // Cambridge → Paris (the canonical haversine worked example,
        // ≈ 404.3 km, initial bearing ≈ 156.2°), with a query point
        // roughly between them.
        Self {
            lat_a: 52.205,
            lon_a: 0.119,
            lat_b: 48.857,
            lon_b: 2.351,
            lat_p: 50.5,
            lon_p: 1.2,
            result: String::new(),
            error: None,
        }
    }
}

/// Draw the Geomatics Workbench right-side panel. A no-op when the
/// `show_geomatics_workbench` toggle is off.
pub fn draw_geomatics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_geomatics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_geomatics_workbench",
        "Geomatics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native geodesic calculations · valenx-geomatics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.geomatics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Point A (deg)").strong());
                    ui.horizontal(|ui| {
                        ui.label("lat");
                        ui.add(egui::DragValue::new(&mut s.lat_a).speed(0.01));
                        ui.label("lon");
                        ui.add(egui::DragValue::new(&mut s.lon_a).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Point B (deg)").strong());
                    ui.horizontal(|ui| {
                        ui.label("lat");
                        ui.add(egui::DragValue::new(&mut s.lat_b).speed(0.01));
                        ui.label("lon");
                        ui.add(egui::DragValue::new(&mut s.lon_b).speed(0.01));
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Query point P (deg)").strong());
                    ui.horizontal(|ui| {
                        ui.label("lat");
                        ui.add(egui::DragValue::new(&mut s.lat_p).speed(0.01));
                        ui.label("lon");
                        ui.add(egui::DragValue::new(&mut s.lon_p).speed(0.01));
                    });
                    ui.label(
                        egui::RichText::new(
                            "cross / along-track measure P against the path A\u{2192}B",
                        )
                        .weak()
                        .small(),
                    );

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("\u{25B6} Compute").strong())
                        .clicked()
                    {
                        run_geomatics(s);
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Geodesics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_geomatics_workbench = false;
    }
}

/// Validate the three points, run the geodesic calculations, and format
/// the readout. Extracted from the draw closure so it is unit-testable.
fn run_geomatics(s: &mut GeomaticsWorkbenchState) {
    s.error = None;

    for (lat, lon, name) in [
        (s.lat_a, s.lon_a, "A"),
        (s.lat_b, s.lon_b, "B"),
        (s.lat_p, s.lon_p, "P"),
    ] {
        if !(lat.is_finite() && (-90.0..=90.0).contains(&lat)) {
            s.error = Some(format!("point {name} latitude must be in [-90, 90]"));
            return;
        }
        if !(lon.is_finite() && (-180.0..=180.0).contains(&lon)) {
            s.error = Some(format!("point {name} longitude must be in [-180, 180]"));
            return;
        }
    }

    // Elevation is ignored by every geodesic routine; carry 0.
    let a = LatLon {
        latitude_deg: s.lat_a,
        longitude_deg: s.lon_a,
        elevation_m: 0.0,
    };
    let b = LatLon {
        latitude_deg: s.lat_b,
        longitude_deg: s.lon_b,
        elevation_m: 0.0,
    };
    let p = LatLon {
        latitude_deg: s.lat_p,
        longitude_deg: s.lon_p,
        elevation_m: 0.0,
    };

    let great_circle_km = haversine_distance(a, b) / 1000.0;
    let init_brg = initial_bearing(a, b);
    let fin_brg = final_bearing(a, b);
    let rhumb_km = rhumb_distance(a, b) / 1000.0;
    let cross_track_m = cross_track_distance(p, a, b);
    let along_track_km = along_track_distance(p, a, b) / 1000.0;

    s.result = format!(
        "point A        : {:.4}\u{00B0}, {:.4}\u{00B0}\n\
         point B        : {:.4}\u{00B0}, {:.4}\u{00B0}\n\
         query point P  : {:.4}\u{00B0}, {:.4}\u{00B0}\n\n\
         great-circle   : {:.3} km   (haversine A\u{2192}B)\n\
         initial bearing: {:.1}\u{00B0}   (A\u{2192}B, cw from north)\n\
         final bearing  : {:.1}\u{00B0}   (arrival heading at B)\n\
         rhumb line     : {:.3} km   (constant bearing A\u{2192}B)\n\
         cross-track    : {:.1} m   (P offset from path; +right/\u{2212}left)\n\
         along-track    : {:.3} km   (P projected onto path from A)",
        s.lat_a,
        s.lon_a,
        s.lat_b,
        s.lon_b,
        s.lat_p,
        s.lon_p,
        great_circle_km,
        init_brg,
        fin_brg,
        rhumb_km,
        cross_track_m,
        along_track_km,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_is_idle() {
        let s = GeomaticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn compute_cambridge_to_paris() {
        let mut s = GeomaticsWorkbenchState::default();
        run_geomatics(&mut s);
        assert!(s.error.is_none());
        assert!(!s.result.is_empty());
        // The readout names each geodesic.
        assert!(s.result.contains("great-circle"));
        assert!(s.result.contains("initial bearing"));
        assert!(s.result.contains("rhumb line"));
        assert!(s.result.contains("cross-track"));
        assert!(s.result.contains("along-track"));
        // Cambridge → Paris great-circle is the canonical ≈ 404.3 km
        // haversine example; recompute via the backend to confirm the
        // defaults are the intended worked example.
        let a = LatLon {
            latitude_deg: 52.205,
            longitude_deg: 0.119,
            elevation_m: 0.0,
        };
        let b = LatLon {
            latitude_deg: 48.857,
            longitude_deg: 2.351,
            elevation_m: 0.0,
        };
        assert!((haversine_distance(a, b) / 1000.0 - 404.3).abs() < 1.0);
        // Initial bearing ≈ 156.2° (south-east), normalised to [0, 360).
        let ib = initial_bearing(a, b);
        assert!((150.0..162.0).contains(&ib));
    }

    #[test]
    fn compute_rejects_out_of_range_coords() {
        for bad in [
            GeomaticsWorkbenchState {
                lat_a: 95.0,
                ..Default::default()
            },
            GeomaticsWorkbenchState {
                lon_b: 200.0,
                ..Default::default()
            },
            GeomaticsWorkbenchState {
                lat_p: f64::NAN,
                ..Default::default()
            },
        ] {
            let mut s = bad;
            run_geomatics(&mut s);
            assert!(s.error.is_some());
            assert!(s.result.is_empty());
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;

    /// Render the whole workbench panel once in a headless egui context.
    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_geomatics_workbench(app, ctx);
        });
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_geomatics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_geomatics_workbench = true;
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_a_result_and_error_without_panic() {
        let mut app = ValenxApp::default();
        app.show_geomatics_workbench = true;
        run_geomatics(&mut app.geomatics);
        app.geomatics.error = Some("invalid coordinates".to_string());
        draw_workbench(&mut app);
    }
}
