//! Closed-form ground-truth validation of the curvature engine.
//!
//! Every check here compares the engine's numeric output against a value that
//! is known analytically, so a regression in the AD pipeline, the metric
//! definitions, or the tensor assembly shows up immediately.

use std::f64::consts::FRAC_PI_2;

use valenx_relativity::{
    curvature_at, kerr, reissner_nordstrom, schwarzschild, KerrNewman, Minkowski,
};

const EQUATOR: f64 = FRAC_PI_2;

// ---- flat baseline ---------------------------------------------------------

#[test]
fn minkowski_is_exactly_flat() {
    let c = curvature_at(&Minkowski, [1.0, 2.0, 3.0, 4.0]).unwrap();
    for a in 0..4 {
        for b in 0..4 {
            for cc in 0..4 {
                for d in 0..4 {
                    assert!(
                        c.riemann[a][b][cc][d].abs() < 1e-12,
                        "Minkowski Riemann nonzero at {a}{b}{cc}{d}: {}",
                        c.riemann[a][b][cc][d]
                    );
                }
            }
        }
    }
    assert!(c.ricci_scalar.abs() < 1e-12);
    assert!(c.kretschmann.abs() < 1e-12);
}

// ---- Schwarzschild: vacuum + Kretschmann = 48 M²/r⁶ ------------------------

#[test]
fn schwarzschild_is_vacuum() {
    let bh = schwarzschild(1.0);
    for &r in &[3.0, 5.0, 10.0, 50.0] {
        let c = curvature_at(&bh, [0.0, r, EQUATOR, 0.0]).unwrap();
        for a in 0..4 {
            for b in 0..4 {
                assert!(
                    c.ricci[a][b].abs() < 1e-7,
                    "Schwarzschild Ricci[{a}][{b}] = {} at r={r}",
                    c.ricci[a][b]
                );
            }
        }
        assert!(
            c.ricci_scalar.abs() < 1e-7,
            "R = {} at r={r}",
            c.ricci_scalar
        );
    }
}

#[test]
fn schwarzschild_kretschmann_matches_closed_form() {
    // K = 48 M² / r⁶ for any mass and radius. Radii are taken in units of M so
    // every sample sits outside the event horizon (r = 2M) regardless of mass.
    for &m in &[0.5, 1.0, 2.0] {
        let bh = schwarzschild(m);
        for &r_over_m in &[3.0, 4.0, 6.0, 10.0, 25.0] {
            let r = r_over_m * m;
            let c = curvature_at(&bh, [0.0, r, EQUATOR, 0.0]).unwrap();
            let expected = 48.0 * m * m / r.powi(6);
            let rel = (c.kretschmann - expected).abs() / expected;
            assert!(
                rel < 1e-8,
                "Kretschmann at M={m}, r={r}: got {}, want {expected} (rel {rel:e})",
                c.kretschmann
            );
        }
    }
}

#[test]
fn schwarzschild_christoffel_matches_analytic_equator() {
    // Hand-derived Schwarzschild Christoffels at M=1, r=5, θ=π/2 (f = 1−2M/r).
    let bh = schwarzschild(1.0);
    let r = 5.0;
    let c = curvature_at(&bh, [0.0, r, EQUATOR, 0.0]).unwrap();
    let g = &c.christoffel; // g[a][b][c] = Γ^a_{bc}; 0=t,1=r,2=θ,3=φ
    let approx = |got: f64, want: f64, name: &str| {
        assert!(
            (got - want).abs() < 1e-9,
            "Γ {name}: got {got}, want {want}"
        );
    };
    // Γ^t_{tr} = M / (r(r−2M)) = 1/15
    approx(g[0][0][1], 1.0 / 15.0, "t_tr");
    approx(g[0][1][0], 1.0 / 15.0, "t_rt (symmetry)");
    // Γ^r_{tt} = M(r−2M)/r³ = 3/125
    approx(g[1][0][0], 3.0 / 125.0, "r_tt");
    // Γ^r_{rr} = −M/(r(r−2M)) = −1/15
    approx(g[1][1][1], -1.0 / 15.0, "r_rr");
    // Γ^r_{θθ} = −(r−2M) = −3
    approx(g[1][2][2], -3.0, "r_thth");
    // Γ^r_{φφ} = −(r−2M) sin²θ = −3
    approx(g[1][3][3], -3.0, "r_phph");
    // Γ^θ_{rθ} = 1/r = 0.2
    approx(g[2][1][2], 0.2, "th_rth");
    approx(g[2][2][1], 0.2, "th_thr (symmetry)");
    // Γ^φ_{rφ} = 1/r = 0.2
    approx(g[3][1][3], 0.2, "ph_rph");
}

#[test]
fn schwarzschild_christoffel_matches_analytic_off_equator() {
    // θ-dependent components at M=1, r=5, θ=π/3.
    let bh = schwarzschild(1.0);
    let theta = std::f64::consts::PI / 3.0;
    let (s, co) = (theta.sin(), theta.cos());
    let c = curvature_at(&bh, [0.0, 5.0, theta, 0.0]).unwrap();
    let g = &c.christoffel;
    let approx = |got: f64, want: f64, name: &str| {
        assert!(
            (got - want).abs() < 1e-9,
            "Γ {name}: got {got}, want {want}"
        );
    };
    // Γ^θ_{φφ} = −sinθ cosθ
    approx(g[2][3][3], -s * co, "th_phph");
    // Γ^φ_{θφ} = cotθ
    approx(g[3][2][3], co / s, "ph_thph");
    // Γ^r_{φφ} = −(r−2M) sin²θ = −3 sin²θ
    approx(g[1][3][3], -3.0 * s * s, "r_phph");
    // Γ^r_{θθ} = −(r−2M) = −3 (θ-independent)
    approx(g[1][2][2], -3.0, "r_thth");
}

// ---- Kerr: vacuum solution -------------------------------------------------

#[test]
fn kerr_is_vacuum() {
    // Kerr is a vacuum solution: Ricci = 0 everywhere outside the horizon,
    // including off-axis / off-equator and for nonzero spin.
    let bh = kerr(1.0, 0.6);
    let points = [
        [0.0, 4.0, EQUATOR, 0.0],
        [0.0, 6.0, std::f64::consts::PI / 3.0, 1.0],
        [0.0, 10.0, std::f64::consts::PI / 4.0, 0.0],
    ];
    for p in points {
        let c = curvature_at(&bh, p).unwrap();
        for a in 0..4 {
            for b in 0..4 {
                assert!(
                    c.ricci[a][b].abs() < 1e-6,
                    "Kerr Ricci[{a}][{b}] = {} at {p:?}",
                    c.ricci[a][b]
                );
            }
        }
        assert!(
            c.ricci_scalar.abs() < 1e-6,
            "Kerr R = {} at {p:?}",
            c.ricci_scalar
        );
    }
}

// ---- charged holes: traceless EM source -> R = 0 but Ricci != 0 ------------

#[test]
fn reissner_nordstrom_has_traceless_ricci() {
    // The electromagnetic stress-energy is traceless, so the Ricci *scalar*
    // vanishes — but the Ricci *tensor* does not (charge curves spacetime).
    let bh = reissner_nordstrom(1.0, 0.5);
    let c = curvature_at(&bh, [0.0, 5.0, EQUATOR, 0.0]).unwrap();
    assert!(
        c.ricci_scalar.abs() < 1e-7,
        "RN R should vanish, got {}",
        c.ricci_scalar
    );
    let max_ricci = (0..4)
        .flat_map(|a| (0..4).map(move |b| (a, b)))
        .map(|(a, b)| c.ricci[a][b].abs())
        .fold(0.0_f64, f64::max);
    assert!(
        max_ricci > 1e-3,
        "RN Ricci tensor should be nonzero (charge sources curvature), max = {max_ricci}"
    );
}

#[test]
fn kerr_newman_reduces_to_special_cases() {
    // Kerr–Newman with a=Q=0 must reproduce Schwarzschild curvature exactly.
    let kn = KerrNewman {
        mass: 1.0,
        spin: 0.0,
        charge: 0.0,
    };
    let sch = schwarzschild(1.0);
    let p = [0.0, 7.0, EQUATOR, 0.0];
    let a = curvature_at(&kn, p).unwrap();
    let b = curvature_at(&sch, p).unwrap();
    assert!((a.kretschmann - b.kretschmann).abs() < 1e-12);
    assert!((a.kretschmann - 48.0 / 7.0_f64.powi(6)).abs() / (48.0 / 7.0_f64.powi(6)) < 1e-8);
}

// ---- coordinate-singularity handling ---------------------------------------

#[test]
fn horizon_is_reported_not_nan() {
    // At the Schwarzschild horizon r = 2M the metric (g_rr) blows up; the
    // engine must return an error, never a silent NaN.
    let bh = schwarzschild(1.0);
    let res = curvature_at(&bh, [0.0, 2.0, EQUATOR, 0.0]);
    assert!(
        res.is_err(),
        "expected a coordinate-singularity error at the horizon"
    );
}

// ---- serialization round-trip ----------------------------------------------

#[test]
fn structs_round_trip_through_json() {
    let bh = KerrNewman {
        mass: 1.0,
        spin: 0.3,
        charge: 0.1,
    };
    let json = serde_json::to_string(&bh).unwrap();
    let back: KerrNewman = serde_json::from_str(&json).unwrap();
    assert_eq!(bh, back);

    let c = curvature_at(&schwarzschild(1.0), [0.0, 6.0, EQUATOR, 0.0]).unwrap();
    let json = serde_json::to_string(&c).unwrap();
    let back: valenx_relativity::Curvature = serde_json::from_str(&json).unwrap();
    assert!((c.kretschmann - back.kretschmann).abs() < 1e-15);
}
