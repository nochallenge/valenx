//! Ready-made example vehicles and run configurations.
//!
//! These are illustrative, not flight data for any real vehicle — they
//! are sized to behave like a generic medium-lift two-stage launcher so
//! the simulator can be exercised end-to-end out of the box.

use crate::config::{AscentConfig, GuidanceMode};
use crate::guidance::GuidanceProgram;
use crate::vehicle::{DragModel, Stage, Vehicle};
use crate::wind::WindModel;

/// A generic two-stage medium-lift launch vehicle (kerolox-class
/// numbers, ~10 t to LEO). Illustrative, not a real vehicle's data.
pub fn two_stage_medium_lift() -> Vehicle {
    Vehicle {
        stages: vec![
            Stage {
                name: "first stage".into(),
                dry_mass: 25_000.0,
                propellant_mass: 410_000.0,
                thrust_vac: 8_200_000.0,
                thrust_sl: 7_600_000.0,
                isp_vac: 311.0,
                isp_sl: 283.0,
            },
            Stage {
                name: "second stage".into(),
                dry_mass: 4_000.0,
                propellant_mass: 100_000.0,
                thrust_vac: 980_000.0,
                thrust_sl: 980_000.0,
                isp_vac: 348.0,
                isp_sl: 348.0,
            },
        ],
        payload_mass: 10_000.0,
        // 3.7 m diameter -> ~10.75 m² frontal area.
        reference_area: 10.75,
        drag: DragModel::generic_launch_vehicle(),
    }
}

/// A launch configuration tuned to drive [`two_stage_medium_lift`] into
/// a bound orbit via a gravity turn.
///
/// With this overpowered vehicle burning all propellant prograde, the
/// open-loop gravity turn reaches a highly *eccentric* orbit (periapsis
/// ~300 km, apoapsis several thousand km). Trimming it to a near-
/// circular LEO requires a coast-to-apoapsis circularisation burn or
/// closed-loop guidance — both documented as future work in the crate
/// docs. The pitch kick sits in the middle of the stable basin
/// (`pk ∈ [11°, 14°]` all reach orbit) so the result is not fragile.
pub fn leo_ascent_config() -> AscentConfig {
    AscentConfig {
        launch_altitude_m: 0.0,
        guidance: GuidanceProgram {
            vertical_rise_time: 20.0,
            pitch_kick_deg: 12.0,
            kick_duration: 5.0,
        },
        time_step: 0.1,
        max_time: 1_500.0,
        sample_interval: 2.0,
        mode: GuidanceMode::OpenLoopGravityTurn,
        wind: WindModel::None,
    }
}

/// A closed-loop insertion configuration that drives
/// [`two_stage_medium_lift`] into a near-circular ~300 km low orbit:
/// ascend → coast to apoapsis → circularise. Uses a gentler pitch kick
/// than [`leo_ascent_config`] so the powered ascent reaches the target
/// apoapsis without lofting, leaving propellant for the circularisation
/// burn.
pub fn leo_insertion_config() -> AscentConfig {
    AscentConfig {
        launch_altitude_m: 0.0,
        guidance: GuidanceProgram {
            vertical_rise_time: 20.0,
            // Centre of the stable basin (pk ∈ [12.6°, 13.2°] all
            // circularise) — a flatter ascent than the open-loop preset
            // so the vehicle arrives at apoapsis with enough horizontal
            // velocity to keep the circularisation burn cheap.
            pitch_kick_deg: 12.9,
            kick_duration: 5.0,
        },
        time_step: 0.1,
        // Allow time for the coast to apoapsis plus the circularisation
        // burn.
        max_time: 3_000.0,
        sample_interval: 2.0,
        mode: GuidanceMode::ClosedLoopInsertion {
            target_altitude_m: 300_000.0,
        },
        wind: WindModel::None,
    }
}
