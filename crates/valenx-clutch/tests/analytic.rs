//! Ground-truth analytic tests for the friction-clutch torque / power
//! models.
//!
//! Every float comparison uses an absolute-difference tolerance, never
//! `==`. The reference numbers are computed independently from the
//! closed-form equations (by hand or from first principles) so the tests
//! validate the implementation against the physics, not against itself.

use valenx_clutch::clutch::rpm_to_rad_per_s;
use valenx_clutch::{ClutchError, ClutchGeometry, FrictionClutch, PressureModel};

/// Absolute tolerance for newton-metre / watt comparisons.
const EPS: f64 = 1e-9;

/// Build a canonical single-plate clutch for reuse across tests.
///
/// 50 mm inner / 100 mm outer face, mu = 0.30, N = 2 surfaces.
fn canonical() -> FrictionClutch {
    let geom = ClutchGeometry::new(50.0, 100.0).expect("valid geometry");
    FrictionClutch::new(geom, 0.30, 2).expect("valid clutch")
}

#[test]
fn uniform_wear_matches_hand_calc() {
    // ri = 0.050 m, ro = 0.100 m, mu = 0.30, N = 2, F = 1000 N.
    // mean radius = (0.100 + 0.050)/2 = 0.075 m.
    // T = mu*F*N*r = 0.30 * 1000 * 2 * 0.075 = 45.0 N*m.
    let clutch = canonical();
    let t = clutch.torque_uniform_wear(1000.0).unwrap();
    assert!((t - 45.0).abs() < EPS, "got {t}");
}

#[test]
fn uniform_pressure_matches_hand_calc() {
    // r_eff = (2/3)(ro^3 - ri^3)/(ro^2 - ri^2)
    //       = (2/3)(1e-3 - 0.125e-3)/(0.01 - 0.0025)
    //       = (2/3)(0.875e-3)/(7.5e-3)
    //       = (2/3) * 0.11666666666666667
    //       = 0.07777777777777778 m.
    // T = 0.30 * 1000 * 2 * r_eff = 600 * r_eff = 46.666666666666664 N*m.
    let clutch = canonical();
    let t = clutch.torque_uniform_pressure(1000.0).unwrap();
    let expected = 600.0 * (2.0 / 3.0) * (0.875e-3) / (7.5e-3);
    assert!((t - expected).abs() < EPS, "got {t}, expected {expected}");
    assert!((t - 46.666_666_666_666_664).abs() < 1e-9, "got {t}");
}

#[test]
fn uniform_wear_is_less_than_uniform_pressure() {
    // Core textbook result: the worn-in mean radius is smaller than the
    // area-weighted mean radius, so uniform-wear torque < uniform-pressure
    // torque for the same geometry / mu / N / F.
    let clutch = canonical();
    let f = 1234.0;
    let tw = clutch.torque_uniform_wear(f).unwrap();
    let tp = clutch.torque_uniform_pressure(f).unwrap();
    assert!(tp > tw, "uw {tw} should be < up {tp}");
    // And strictly so by a finite margin (not a floating tie).
    assert!((tp - tw) > 1e-3, "margin {} too small", tp - tw);
}

#[test]
fn torque_scales_linearly_with_friction_coefficient() {
    // Doubling mu doubles the torque (everything else fixed).
    let geom = ClutchGeometry::new(50.0, 100.0).unwrap();
    let c1 = FrictionClutch::new(geom, 0.20, 2).unwrap();
    let c2 = FrictionClutch::new(geom, 0.40, 2).unwrap();
    let t1 = c1.torque_uniform_wear(1000.0).unwrap();
    let t2 = c2.torque_uniform_wear(1000.0).unwrap();
    assert!((t2 - 2.0 * t1).abs() < EPS, "t1 {t1}, t2 {t2}");
}

#[test]
fn torque_scales_linearly_with_clamp_force() {
    // Tripling F triples the torque (both models).
    let clutch = canonical();
    for model in [PressureModel::UniformWear, PressureModel::UniformPressure] {
        let t1 = clutch.torque(model, 500.0).unwrap();
        let t3 = clutch.torque(model, 1500.0).unwrap();
        assert!((t3 - 3.0 * t1).abs() < EPS, "{model:?}: t1 {t1}, t3 {t3}");
    }
}

#[test]
fn torque_scales_linearly_with_surface_count() {
    // More friction surfaces -> proportionally more torque.
    // N = 4 should give exactly twice the torque of N = 2.
    let geom = ClutchGeometry::new(50.0, 100.0).unwrap();
    let single = FrictionClutch::new(geom, 0.30, 2).unwrap();
    let multi = FrictionClutch::new(geom, 0.30, 4).unwrap();
    let ts = single.torque_uniform_wear(1000.0).unwrap();
    let tm = multi.torque_uniform_wear(1000.0).unwrap();
    assert!((tm - 2.0 * ts).abs() < EPS, "single {ts}, multi {tm}");
    // Strictly monotone in N.
    assert!(tm > ts, "more surfaces must give more torque");
}

#[test]
fn more_surfaces_strictly_increase_torque_monotonically() {
    let geom = ClutchGeometry::new(40.0, 90.0).unwrap();
    let mut last = 0.0;
    for n in 1..=8u32 {
        let clutch = FrictionClutch::new(geom, 0.25, n).unwrap();
        let t = clutch.torque_uniform_pressure(800.0).unwrap();
        if n > 1 {
            assert!(
                t > last,
                "N={n}: torque {t} not greater than previous {last}"
            );
        }
        // N appears linearly: T(n) = n * T(1).
        let t1 = FrictionClutch::new(geom, 0.25, 1)
            .unwrap()
            .torque_uniform_pressure(800.0)
            .unwrap();
        assert!(
            (t - n as f64 * t1).abs() < EPS,
            "N={n}: {t} != {} ",
            n as f64 * t1
        );
        last = t;
    }
}

#[test]
fn torque_increases_with_radius() {
    // A larger outer radius (bigger lever arm) gives more torque for the
    // same clamp force, mu, and N.
    let small = FrictionClutch::new(ClutchGeometry::new(50.0, 80.0).unwrap(), 0.30, 2).unwrap();
    let large = FrictionClutch::new(ClutchGeometry::new(50.0, 120.0).unwrap(), 0.30, 2).unwrap();
    let ts = small.torque_uniform_wear(1000.0).unwrap();
    let tl = large.torque_uniform_wear(1000.0).unwrap();
    assert!(tl > ts, "larger face {tl} should beat smaller {ts}");
}

#[test]
fn uniform_wear_mean_radius_is_arithmetic_mean() {
    let geom = ClutchGeometry::new(50.0, 100.0).unwrap();
    // (0.050 + 0.100)/2 = 0.075 m.
    assert!((geom.mean_radius_uniform_wear_m() - 0.075).abs() < EPS);
}

#[test]
fn uniform_pressure_mean_radius_exceeds_uniform_wear_mean_radius() {
    // Geometric fact: the centroidal mean radius is strictly larger than
    // the arithmetic mean radius for any non-degenerate annulus.
    for (ri, ro) in [(10.0, 20.0), (50.0, 100.0), (5.0, 200.0), (99.0, 100.0)] {
        let g = ClutchGeometry::new(ri, ro).unwrap();
        let rw = g.mean_radius_uniform_wear_m();
        let rp = g.mean_radius_uniform_pressure_m();
        assert!(rp > rw, "ri={ri} ro={ro}: rp {rp} not > rw {rw}");
        // Both lie inside the annulus.
        assert!(rw > g.inner_radius_m() && rw < g.outer_radius_m());
        assert!(rp > g.inner_radius_m() && rp < g.outer_radius_m());
    }
}

#[test]
fn power_is_torque_times_omega() {
    // P = T * omega, verified against an independent torque evaluation.
    let clutch = canonical();
    let f = 2000.0;
    let omega = 150.0; // rad/s
    let t = clutch.torque(PressureModel::UniformWear, f).unwrap();
    let p = clutch.power(PressureModel::UniformWear, f, omega).unwrap();
    assert!((p - t * omega).abs() < EPS, "t {t}, omega {omega}, p {p}");
}

#[test]
fn power_scales_linearly_with_omega() {
    let clutch = canonical();
    let f = 2000.0;
    let p1 = clutch
        .power(PressureModel::UniformPressure, f, 100.0)
        .unwrap();
    let p2 = clutch
        .power(PressureModel::UniformPressure, f, 200.0)
        .unwrap();
    assert!((p2 - 2.0 * p1).abs() < EPS, "p1 {p1}, p2 {p2}");
}

#[test]
fn power_at_zero_speed_is_zero() {
    let clutch = canonical();
    let p = clutch
        .power(PressureModel::UniformWear, 1000.0, 0.0)
        .unwrap();
    assert!(p.abs() < EPS, "static clutch transmits no power, got {p}");
}

#[test]
fn dispatch_torque_matches_named_methods() {
    let clutch = canonical();
    let f = 777.0;
    let tw = clutch.torque_uniform_wear(f).unwrap();
    let tp = clutch.torque_uniform_pressure(f).unwrap();
    assert!((clutch.torque(PressureModel::UniformWear, f).unwrap() - tw).abs() < EPS);
    assert!((clutch.torque(PressureModel::UniformPressure, f).unwrap() - tp).abs() < EPS);
}

#[test]
fn rpm_conversion_matches_definition() {
    // 60 rpm = 1 rev/s = 2*pi rad/s.
    let w = rpm_to_rad_per_s(60.0).unwrap();
    assert!((w - std::f64::consts::TAU).abs() < EPS, "got {w}");
    // 3000 rpm = 50 rev/s = 100*pi rad/s.
    let w2 = rpm_to_rad_per_s(3000.0).unwrap();
    assert!((w2 - 100.0 * std::f64::consts::PI).abs() < 1e-9, "got {w2}");
}

#[test]
fn torque_curve_is_linear_through_origin() {
    let clutch = canonical();
    let forces = [0.0, 250.0, 500.0, 750.0, 1000.0];
    let curve = clutch
        .torque_curve(PressureModel::UniformWear, &forces)
        .unwrap();
    assert_eq!(curve.len(), forces.len());
    // Slope mu*N*r_eff; verify each point equals slope * force.
    let slope = clutch.torque_uniform_wear(1.0).unwrap();
    for (i, &f) in forces.iter().enumerate() {
        assert!(
            (curve[i] - slope * f).abs() < EPS,
            "i={i}: {} vs {}",
            curve[i],
            slope * f
        );
    }
    // Zero clamp force -> zero torque.
    assert!(curve[0].abs() < EPS);
}

#[test]
fn zero_clamp_force_gives_zero_torque() {
    let clutch = canonical();
    assert!(clutch.torque_uniform_wear(0.0).unwrap().abs() < EPS);
    assert!(clutch.torque_uniform_pressure(0.0).unwrap().abs() < EPS);
}

// --- validation of the error paths --------------------------------------

#[test]
fn rejects_inverted_radii() {
    let err = ClutchGeometry::new(100.0, 50.0).unwrap_err();
    assert!(matches!(err, ClutchError::InvertedRadii { .. }));
    assert_eq!(err.code(), "clutch.inverted-radii");
}

#[test]
fn rejects_equal_radii() {
    // ri == ro is degenerate (empty annulus): also InvertedRadii.
    let err = ClutchGeometry::new(75.0, 75.0).unwrap_err();
    assert!(matches!(err, ClutchError::InvertedRadii { .. }));
}

#[test]
fn rejects_non_positive_radius() {
    assert!(matches!(
        ClutchGeometry::new(0.0, 100.0).unwrap_err(),
        ClutchError::InvalidParameter {
            name: "inner_mm",
            ..
        }
    ));
    assert!(matches!(
        ClutchGeometry::new(-5.0, 100.0).unwrap_err(),
        ClutchError::InvalidParameter {
            name: "inner_mm",
            ..
        }
    ));
}

#[test]
fn rejects_non_finite_radius() {
    assert!(matches!(
        ClutchGeometry::new(f64::NAN, 100.0).unwrap_err(),
        ClutchError::InvalidParameter { .. }
    ));
    assert!(matches!(
        ClutchGeometry::new(50.0, f64::INFINITY).unwrap_err(),
        ClutchError::InvalidParameter { .. }
    ));
}

#[test]
fn rejects_non_positive_friction_coefficient() {
    let geom = ClutchGeometry::new(50.0, 100.0).unwrap();
    assert!(matches!(
        FrictionClutch::new(geom, 0.0, 2).unwrap_err(),
        ClutchError::InvalidParameter { name: "mu", .. }
    ));
    assert!(matches!(
        FrictionClutch::new(geom, -0.3, 2).unwrap_err(),
        ClutchError::InvalidParameter { name: "mu", .. }
    ));
}

#[test]
fn rejects_zero_surface_count() {
    let geom = ClutchGeometry::new(50.0, 100.0).unwrap();
    let err = FrictionClutch::new(geom, 0.30, 0).unwrap_err();
    assert!(matches!(err, ClutchError::InvalidSurfaceCount(_)));
    assert_eq!(err.code(), "clutch.invalid-surface-count");
}

#[test]
fn rejects_negative_clamp_force() {
    let clutch = canonical();
    assert!(matches!(
        clutch.torque_uniform_wear(-1.0).unwrap_err(),
        ClutchError::InvalidParameter {
            name: "clamp_force_n",
            ..
        }
    ));
}

#[test]
fn rejects_negative_speed() {
    let clutch = canonical();
    assert!(matches!(
        clutch
            .power(PressureModel::UniformWear, 1000.0, -10.0)
            .unwrap_err(),
        ClutchError::InvalidParameter {
            name: "omega_rad_per_s",
            ..
        }
    ));
    assert!(matches!(
        rpm_to_rad_per_s(-1.0).unwrap_err(),
        ClutchError::InvalidParameter { name: "rpm", .. }
    ));
}

#[test]
fn from_metres_and_mm_constructors_agree() {
    let a = ClutchGeometry::new(50.0, 100.0).unwrap();
    let b = ClutchGeometry::from_metres(0.050, 0.100).unwrap();
    assert!((a.inner_radius_m() - b.inner_radius_m()).abs() < EPS);
    assert!((a.outer_radius_m() - b.outer_radius_m()).abs() < EPS);
}

#[test]
fn serde_round_trip_preserves_clutch() {
    let clutch = canonical();
    let json = serde_json::to_string(&clutch).unwrap();
    let back: FrictionClutch = serde_json::from_str(&json).unwrap();
    assert_eq!(clutch, back);
    // Behaviour survives the round trip too.
    let f = 1500.0;
    assert!(
        (clutch.torque_uniform_wear(f).unwrap() - back.torque_uniform_wear(f).unwrap()).abs() < EPS
    );
}

#[test]
fn worked_textbook_example_single_plate() {
    // Worked example: single-plate dry clutch, ri = 75 mm, ro = 150 mm,
    // mu = 0.25, clamp force F = 4 kN, N = 2 faces, at 2000 rpm.
    //
    // Uniform wear:
    //   r_mean = (0.150 + 0.075)/2 = 0.1125 m
    //   T = 0.25 * 4000 * 2 * 0.1125 = 225.0 N*m
    //   omega = 2*pi*2000/60 = 209.4395102393... rad/s
    //   P = T*omega = 225 * 209.4395102393... = 47123.889803... W
    let geom = ClutchGeometry::new(75.0, 150.0).unwrap();
    let clutch = FrictionClutch::new(geom, 0.25, 2).unwrap();
    let t = clutch.torque_uniform_wear(4000.0).unwrap();
    assert!((t - 225.0).abs() < EPS, "T = {t}");

    let omega = rpm_to_rad_per_s(2000.0).unwrap();
    let expected_omega = 2.0 * std::f64::consts::PI * 2000.0 / 60.0;
    assert!((omega - expected_omega).abs() < EPS);

    let p = clutch
        .power(PressureModel::UniformWear, 4000.0, omega)
        .unwrap();
    assert!((p - 225.0 * expected_omega).abs() < 1e-6, "P = {p}");
    // Sanity on the magnitude: ~47.1 kW.
    assert!((p - 47_123.889_803_846_9).abs() < 1e-3, "P = {p}");
}
