//! Closed-form ground-truth validation of observables, thermodynamics and the
//! SI unit conversions.

use std::f64::consts::{FRAC_PI_2, PI};

use valenx_relativity::observables::OrbitSense::{Prograde, Retrograde};
use valenx_relativity::{
    ergosphere_radius, gravitational_redshift, horizons, isco, kerr, photon_sphere,
    reissner_nordstrom, schwarzschild, shadow_radius, thermodynamics, units, KerrNewman,
};

fn close(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() <= tol,
        "{what}: got {a}, want {b} (tol {tol})"
    );
}
fn rel(a: f64, b: f64, tol: f64, what: &str) {
    assert!(
        (a - b).abs() / b.abs() <= tol,
        "{what}: got {a}, want {b} (rel {})",
        (a - b).abs() / b.abs()
    );
}

// ---- horizons --------------------------------------------------------------

#[test]
fn horizon_radii() {
    let h = horizons(&schwarzschild(1.0)).unwrap();
    close(h.outer, 2.0, 1e-12, "Schwarzschild r+");
    close(h.inner, 0.0, 1e-12, "Schwarzschild r-");

    // Kerr a=0.6: r± = 1 ± √(1−0.36) = 1 ± 0.8
    let h = horizons(&kerr(1.0, 0.6)).unwrap();
    close(h.outer, 1.8, 1e-12, "Kerr r+");
    close(h.inner, 0.2, 1e-12, "Kerr r-");

    // Reissner–Nordström Q=0.6: same discriminant as Kerr a=0.6.
    let h = horizons(&reissner_nordstrom(1.0, 0.6)).unwrap();
    close(h.outer, 1.8, 1e-12, "RN r+");

    // Extremal Kerr a=M: r+ = r- = M.
    let h = horizons(&kerr(1.0, 1.0)).unwrap();
    close(h.outer, 1.0, 1e-12, "extremal Kerr r+");
    close(h.inner, 1.0, 1e-12, "extremal Kerr r-");

    // Super-extremal: no horizon -> error.
    assert!(horizons(&kerr(1.0, 1.5)).is_err());
    assert!(horizons(&KerrNewman {
        mass: 1.0,
        spin: 0.8,
        charge: 0.8
    })
    .is_err());
    assert!(horizons(&schwarzschild(-1.0)).is_err());
}

#[test]
fn ergosphere_touches_horizon_at_poles() {
    let bh = kerr(1.0, 0.6);
    // Equator: r_E = M + √(M²) = 2M (Kerr).
    close(
        ergosphere_radius(&bh, FRAC_PI_2).unwrap(),
        2.0,
        1e-12,
        "Kerr ergo equator",
    );
    // Pole: r_E = M + √(M² − a²) = r+ .
    let rplus = horizons(&bh).unwrap().outer;
    close(
        ergosphere_radius(&bh, 0.0).unwrap(),
        rplus,
        1e-12,
        "Kerr ergo pole = r+",
    );
}

// ---- photon sphere ---------------------------------------------------------

#[test]
fn photon_sphere_values() {
    close(
        photon_sphere(&schwarzschild(1.0), Prograde).unwrap(),
        3.0,
        1e-12,
        "Schwarzschild photon",
    );

    // Reissner–Nordström: (3M + √(9M² − 8Q²))/2.
    let q: f64 = 0.5;
    let want = (3.0 + (9.0 - 8.0 * q * q).sqrt()) / 2.0;
    close(
        photon_sphere(&reissner_nordstrom(1.0, q), Prograde).unwrap(),
        want,
        1e-12,
        "RN photon",
    );

    // Extremal Kerr: prograde photon orbit at r=M, retrograde at r=4M.
    close(
        photon_sphere(&kerr(1.0, 1.0), Prograde).unwrap(),
        1.0,
        1e-9,
        "extremal Kerr photon prograde",
    );
    close(
        photon_sphere(&kerr(1.0, 1.0), Retrograde).unwrap(),
        4.0,
        1e-9,
        "extremal Kerr photon retrograde",
    );

    // Spin AND charge -> no closed form.
    assert!(photon_sphere(
        &KerrNewman {
            mass: 1.0,
            spin: 0.3,
            charge: 0.3
        },
        Prograde
    )
    .is_err());
}

// ---- ISCO ------------------------------------------------------------------

#[test]
fn isco_values() {
    close(
        isco(&schwarzschild(1.0), Prograde).unwrap(),
        6.0,
        1e-9,
        "Schwarzschild ISCO",
    );
    // Extremal Kerr: prograde ISCO at r=M, retrograde at r=9M.
    close(
        isco(&kerr(1.0, 1.0), Prograde).unwrap(),
        1.0,
        1e-9,
        "extremal Kerr ISCO prograde",
    );
    close(
        isco(&kerr(1.0, 1.0), Retrograde).unwrap(),
        9.0,
        1e-9,
        "extremal Kerr ISCO retrograde",
    );
    // Prograde ISCO shrinks with spin, retrograde grows.
    let a = 0.9;
    assert!(isco(&kerr(1.0, a), Prograde).unwrap() < 6.0);
    assert!(isco(&kerr(1.0, a), Retrograde).unwrap() > 6.0);
    // Charge -> no closed form.
    assert!(isco(&reissner_nordstrom(1.0, 0.5), Prograde).is_err());
}

// ---- shadow ----------------------------------------------------------------

#[test]
fn shadow_radius_schwarzschild_is_sqrt27() {
    let b = shadow_radius(&schwarzschild(1.0)).unwrap();
    close(b, 27.0_f64.sqrt(), 1e-9, "Schwarzschild shadow = √27 M");
    close(b, 3.0 * 3.0_f64.sqrt(), 1e-9, "= 3√3 M");
    // Rotating -> not a single radius.
    assert!(shadow_radius(&kerr(1.0, 0.5)).is_err());
}

// ---- redshift --------------------------------------------------------------

#[test]
fn redshift_to_infinity() {
    // Emitted at r=3M, observed far away: 1+z = 1/√(1−2/3) = √3.
    let z = gravitational_redshift(&schwarzschild(1.0), 3.0, 1.0e9).unwrap();
    rel(z, 3.0_f64.sqrt(), 1e-6, "Schwarzschild redshift r=3M -> ∞");
    // No static observer at/inside the horizon.
    assert!(gravitational_redshift(&schwarzschild(1.0), 2.0, 10.0).is_err());
}

// ---- thermodynamics --------------------------------------------------------

#[test]
fn schwarzschild_thermodynamics() {
    let t = thermodynamics(&schwarzschild(1.0)).unwrap();
    close(t.surface_gravity, 0.25, 1e-12, "κ = 1/4M");
    close(t.hawking_temperature, 1.0 / (8.0 * PI), 1e-12, "T = 1/8πM");
    close(t.horizon_area, 16.0 * PI, 1e-12, "A = 16πM²");
    close(t.entropy, 4.0 * PI, 1e-12, "S = 4πM²");
    close(t.horizon_angular_velocity, 0.0, 1e-12, "Ω = 0");
}

#[test]
fn temperature_scales_inverse_mass() {
    let t1 = thermodynamics(&schwarzschild(1.0))
        .unwrap()
        .hawking_temperature;
    let t2 = thermodynamics(&schwarzschild(2.0))
        .unwrap()
        .hawking_temperature;
    rel(t1 / t2, 2.0, 1e-12, "T ∝ 1/M");
}

#[test]
fn extremal_kerr_has_zero_temperature() {
    let t = thermodynamics(&kerr(1.0, 1.0)).unwrap();
    close(t.surface_gravity, 0.0, 1e-12, "extremal κ = 0");
    close(t.hawking_temperature, 0.0, 1e-12, "extremal T = 0");
}

// ---- SI units --------------------------------------------------------------

#[test]
fn si_conversions() {
    // The Sun's Schwarzschild radius is ≈ 2.95 km.
    rel(
        units::schwarzschild_radius_km(1.0),
        2.953,
        2e-3,
        "Sun r_s km",
    );
    // A solar-mass hole's Hawking temperature ≈ 6.17e-8 K.
    rel(
        units::hawking_temperature_kelvin(1.0),
        6.17e-8,
        5e-3,
        "Sun Hawking T",
    );
    // Evaporation time ≈ 2.1e67 years.
    rel(
        units::evaporation_time_years(1.0),
        2.1e67,
        5e-2,
        "Sun evaporation time",
    );
    // Hawking temperature scales as 1/M.
    rel(
        units::hawking_temperature_kelvin(1.0) / units::hawking_temperature_kelvin(10.0),
        10.0,
        1e-9,
        "T ∝ 1/M (SI)",
    );
}

// ---- serialization ---------------------------------------------------------

#[test]
fn thermo_round_trips() {
    let t = thermodynamics(&kerr(1.0, 0.5)).unwrap();
    let j = serde_json::to_string(&t).unwrap();
    let back: valenx_relativity::Thermodynamics = serde_json::from_str(&j).unwrap();
    assert_eq!(t, back);
}
