//! Ground-truth validation of the geodesic integrator: conservation laws,
//! circular orbits, light deflection, perihelion precession and photon capture.

use valenx_relativity::geodesics::{Kind, StopReason};
use valenx_relativity::{
    angular_momentum, deflection_angle, energy, equatorial_state, integrate_geodesic,
    light_deflection_weak_field, norm, orbit_precession, perihelion_advance_per_orbit,
    schwarzschild, GeodesicOptions,
};

/// A null geodesic must conserve energy, angular momentum and its (zero) norm.
#[test]
fn null_geodesic_conserves_invariants() {
    let bh = schwarzschild(1.0);
    let init = equatorial_state(&bh, 50.0, 1.0, 10.0, Kind::Null, true).unwrap();
    let e0 = energy(&bh, &init);
    let l0 = angular_momentum(&bh, &init);
    let opts = GeodesicOptions {
        r_capture: 2.1,
        r_escape: 300.0,
        tol: 1e-12,
        step: 1.0,
        ..Default::default()
    };
    let traj = integrate_geodesic(&bh, init, opts).unwrap();
    assert!(traj.states.len() > 10);
    for st in &traj.states {
        assert!((energy(&bh, st) - e0).abs() < 1e-6, "E drift");
        assert!((angular_momentum(&bh, st) - l0).abs() < 1e-6, "L drift");
        assert!(
            norm(&bh, st).abs() < 1e-6,
            "null norm drift: {}",
            norm(&bh, st)
        );
    }
}

/// A circular timelike orbit at r=10M (with the textbook E, L) stays at r=10M
/// and keeps its timelike norm of −1.
#[test]
fn circular_orbit_stays_circular() {
    let bh = schwarzschild(1.0);
    let r: f64 = 10.0;
    // Schwarzschild circular orbit: L = √(M r²/(r−3M)), E = (1−2M/r)/√(1−3M/r).
    let l = (r * r / (r - 3.0)).sqrt();
    let e = (1.0 - 2.0 / r) / (1.0 - 3.0 / r).sqrt();
    let init = equatorial_state(&bh, r, e, l, Kind::Timelike, false).unwrap();
    let opts = GeodesicOptions {
        r_capture: 2.1,
        r_escape: 1e9,
        max_lambda: 200.0,
        tol: 1e-12,
        step: 0.5,
        ..Default::default()
    };
    let traj = integrate_geodesic(&bh, init, opts).unwrap();
    for st in &traj.states {
        assert!((st.x[1] - r).abs() < 1e-3, "r drifted to {}", st.x[1]);
        assert!((norm(&bh, st) + 1.0).abs() < 1e-6, "timelike norm");
    }
    assert!(traj.last().x[3] > 6.0, "should sweep ~a full revolution");
}

/// Numerically integrated light deflection approaches 4M/b in the weak field.
#[test]
fn light_deflection_matches_weak_field() {
    let bh = schwarzschild(1.0);
    let b = 1000.0;
    let got = deflection_angle(&bh, b).unwrap();
    let want = light_deflection_weak_field(1.0, b); // 0.004 rad
    let rel = (got - want).abs() / want;
    assert!(rel < 0.02, "deflection got {got}, want {want} (rel {rel})");
}

/// Numerically integrated perihelion precession approaches 6πM/p in the weak
/// field. Mercury's 43″/century is this same formula scaled to its orbit.
#[test]
fn perihelion_precession_matches_weak_field() {
    let bh = schwarzschild(1.0);
    let (r1, r2) = (666.666_67, 2000.0); // p = 2 r1 r2/(r1+r2) = 1000 M
    let got = orbit_precession(&bh, r1, r2).unwrap();
    let p = 2.0 * r1 * r2 / (r1 + r2);
    let want = perihelion_advance_per_orbit(1.0, p); // 6π/1000 ≈ 0.018850
    let rel = (got - want).abs() / want;
    assert!(rel < 0.03, "precession got {got}, want {want} (rel {rel})");
}

/// Photons below the critical impact parameter b_crit = 3√3 M ≈ 5.196 are
/// captured; those above escape.
#[test]
fn photon_capture_threshold() {
    let bh = schwarzschild(1.0);
    let opts = |esc: f64| GeodesicOptions {
        r_capture: 2.05,
        r_escape: esc,
        tol: 1e-11,
        step: 1.0,
        ..Default::default()
    };
    // b = 4 < b_crit -> captured.
    let captured = equatorial_state(&bh, 100.0, 1.0, 4.0, Kind::Null, true).unwrap();
    let traj = integrate_geodesic(&bh, captured, opts(200.0)).unwrap();
    assert_eq!(traj.stop, StopReason::Captured, "b=4 should be captured");
    // b = 7 > b_crit -> escapes.
    let escapes = equatorial_state(&bh, 100.0, 1.0, 7.0, Kind::Null, true).unwrap();
    let traj = integrate_geodesic(&bh, escapes, opts(200.0)).unwrap();
    assert_eq!(traj.stop, StopReason::Escaped, "b=7 should escape");
}
