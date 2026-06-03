//! The Astro / Launch workbench's Run action.
//!
//! Unlike the CFD-side wind tunnel — whose steady RANS solve runs for
//! seconds-to-minutes on a background thread — the launch ascent is a
//! **bounded fixed-step RK4 integration** with a hard step cap baked
//! into [`valenx_astro::config::AscentConfig::validate`]. A LEO ascent
//! is a few thousand steps and completes in well under a frame's worth
//! of time, so the workbench runs it **synchronously on click**: build
//! the vehicle + config from the form, call
//! [`valenx_astro::simulate_ascent`], and store the result (or surface
//! the error). No threads, no channels, no `pump_*` — the panel reads
//! `app.astro.last_result` the same frame.
//!
//! This keeps the v1 simple and fully testable: the run path is a plain
//! `&mut ValenxApp -> ()` function a `#[test]` drives directly.

use crate::ValenxApp;
use valenx_astro::simulate_ascent;

use super::model;

/// Validate the form, run the ascent synchronously, and store the
/// result. On error the message lands in `app.astro.error` (shown in a
/// neutral red status line) and `last_result` is left untouched-ish —
/// it is cleared so the stale result of a previous run never reads as
/// the outcome of this failed one.
///
/// Records the form on the undo stack first, so a later Ctrl+Z rewinds
/// the user to the settings of the last run (matches the aero
/// workbench's `record_form` on Run).
pub fn run_ascent(app: &mut ValenxApp) {
    app.astro.error = None;
    app.astro.record_form();

    let vehicle = app.astro.ascent.build_vehicle();
    let config = app.astro.ascent.build_config();

    // simulate_ascent validates the vehicle + config internally and
    // returns a Result — never `.unwrap()` a runtime value.
    match simulate_ascent(&vehicle, &config) {
        Ok(result) => {
            app.astro.status = format!(
                "Ascent complete \u{2014} {:?}, apoapsis {:.0} km, periapsis {:.0} km, \
                 \u{0394}v budget {:.0} m/s",
                result.outcome,
                result.apoapsis_km(),
                result.periapsis_km(),
                result.ideal_delta_v,
            );
            app.astro.last_result = Some(Box::new(result));
        }
        Err(e) => {
            app.astro.status = "Ascent run failed".to_string();
            app.astro.error = Some(model::friendly_error(&e));
            app.astro.last_result = None;
        }
    }
}

/// Public entry-point for the Ctrl+R shortcut handler — mirrors the
/// aero workbench's `start_run_from_shortcut`. The astro run is
/// synchronous, so there is no in-flight guard to honour; this is a
/// thin alias kept for symmetry with the dispatcher.
pub fn run_ascent_from_shortcut(app: &mut ValenxApp) {
    run_ascent(app);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::model::GuidanceChoice;

    #[test]
    fn default_run_reaches_orbit_and_stores_a_result() {
        // The default form is the medium-lift preset — a click must
        // produce a stored result with no error.
        let mut app = ValenxApp::default();
        run_ascent(&mut app);
        assert!(app.astro.error.is_none(), "unexpected error: {:?}", app.astro.error);
        let r = app.astro.last_result.as_ref().expect("a result was stored");
        assert!(r.reached_space, "apoapsis only {:.1} km", r.apoapsis_km());
        assert!(!app.astro.status.is_empty());
    }

    #[test]
    fn closed_loop_run_inserts_near_circular() {
        // The closed-loop insertion mode should circularise near the
        // target altitude — a real LEO, low eccentricity.
        let mut app = ValenxApp::default();
        app.astro.ascent.guidance = GuidanceChoice::ClosedLoopInsertion;
        app.astro.ascent.target_altitude_km = 300.0;
        // The insertion preset uses a gentler pitch kick than the
        // open-loop default; match it so the ascent reaches apoapsis
        // without lofting.
        app.astro.ascent.pitch_kick_deg = 12.9;
        run_ascent(&mut app);
        assert!(app.astro.error.is_none(), "unexpected error: {:?}", app.astro.error);
        let r = app.astro.last_result.as_ref().expect("a result was stored");
        assert!(r.reached_orbit, "periapsis only {:.1} km", r.periapsis_km());
        assert!(r.orbit.eccentricity < 0.05, "ecc {:.4}", r.orbit.eccentricity);
    }

    #[test]
    fn invalid_vehicle_surfaces_error_without_panic() {
        // An empty stage stack is a non-physical vehicle — the run must
        // surface a clean error, clear any stale result, and not panic.
        let mut app = ValenxApp::default();
        // First do a good run so there's a stale result to clear.
        run_ascent(&mut app);
        assert!(app.astro.last_result.is_some());
        // Now break the vehicle and re-run.
        app.astro.ascent.stages.clear();
        run_ascent(&mut app);
        assert!(app.astro.error.is_some(), "empty stack should error");
        assert!(app.astro.last_result.is_none(), "stale result must be cleared");
    }

    #[test]
    fn non_physical_stage_surfaces_error() {
        // A zero-Isp stage is rejected by the backend validation — the
        // UI shows the error, never a NaN result or a panic.
        let mut app = ValenxApp::default();
        app.astro.ascent.stages[0].isp_vac = 0.0;
        run_ascent(&mut app);
        assert!(app.astro.error.is_some(), "zero Isp should error");
        assert!(app.astro.last_result.is_none());
    }

    #[test]
    fn run_records_form_on_the_undo_stack() {
        // Each run snapshots the form so Ctrl+Z rewinds the last run's
        // settings — mirrors the aero workbench.
        let mut app = ValenxApp::default();
        assert!(!app.astro.can_undo());
        run_ascent(&mut app);
        assert!(app.astro.can_undo());
    }
}
