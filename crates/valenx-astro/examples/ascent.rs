//! Fly the bundled two-stage example vehicle to orbit and print a
//! flight report.
//!
//! ```sh
//! cargo run -p valenx-astro --example ascent
//! ```

use valenx_astro::{presets, simulate_ascent};

fn main() {
    let vehicle = presets::two_stage_medium_lift();
    let config = presets::leo_ascent_config();

    let r = simulate_ascent(&vehicle, &config).expect("valid ascent case");

    println!("== Valenx ascent simulation ==");
    println!("liftoff mass      : {:>10.0} kg", r.liftoff_mass);
    println!("ideal Δv budget   : {:>10.0} m/s", r.ideal_delta_v);
    println!("outcome           : {:?}", r.outcome);
    println!();
    println!("apoapsis          : {:>10.0} km", r.apoapsis_km());
    println!("periapsis         : {:>10.0} km", r.periapsis_km());
    println!("eccentricity      : {:>10.3}", r.orbit.eccentricity);
    println!(
        "reached space     : {}   reached orbit: {}",
        r.reached_space, r.reached_orbit
    );
    println!();
    println!(
        "max dynamic press : {:>10.1} kPa @ {:.1} km",
        r.max_dynamic_pressure / 1000.0,
        r.max_q_altitude_m / 1000.0
    );
    println!("max acceleration  : {:>10.1} g", r.max_acceleration_g);
    println!(
        "MECO              : t = {:.0} s, {:.0} m/s, {:.0} km",
        r.final_time,
        r.final_speed_inertial,
        r.final_altitude_m / 1000.0
    );

    println!("\nflight events:");
    for e in &r.events {
        println!(
            "  t = {:>6.1} s   {:>8.1} km   {:>6.0} m/s   {}",
            e.time,
            e.altitude_m / 1000.0,
            e.speed,
            e.kind
        );
    }

    // Closed-loop insertion: ascent -> coast -> circularise into a
    // near-circular LEO.
    let ins =
        simulate_ascent(&vehicle, &presets::leo_insertion_config()).expect("valid insertion case");
    println!("\n== Closed-loop orbital insertion (target 300 km circular) ==");
    println!(
        "apoapsis {:.0} km, periapsis {:.0} km, eccentricity {:.4}  (reached_orbit: {})",
        ins.apoapsis_km(),
        ins.periapsis_km(),
        ins.orbit.eccentricity,
        ins.reached_orbit
    );
    for e in &ins.events {
        println!(
            "  t = {:>6.1} s   {:>8.1} km   {:>6.0} m/s   {}",
            e.time,
            e.altitude_m / 1000.0,
            e.speed,
            e.kind
        );
    }
}
