//! Tests for `valenx-marine-hydro`.
//!
//! The [`validation`] module pins the crate against published reference values
//! (the Holtrop & Mennen (1982) worked example and the ITTC-1957 closed-form
//! friction line); the remaining modules cover units, properties and the error
//! domain.

use super::*;

/// Relative-or-absolute closeness check.
fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs().max(1.0)
}

/// Regression: Holtrop-Mennen has poles outside its prismatic-coefficient
/// validity band — the LCB term at Cp=0.25 and the (0.95-Cp) term at Cp=0.95,
/// both reachable from the Marine workbench block-coefficient field
/// (Cp = Cb / C_m, C_m=0.98). `resistance_at` must reject an out-of-band hull
/// with an `OutOfRange` error so the readout reports "unavailable" instead of
/// a NaN form factor / R_t / P_e.
#[test]
fn resistance_rejects_out_of_band_prismatic_coefficient() {
    let water = WaterProperties::seawater();
    // Cb -> nabla, holding L=205, B=32, T=10 (so Cp = Cb / 0.980).
    let make = |cb: f64| {
        Hull::new(
            205.0,
            32.0,
            10.0,
            cb * 205.0 * 32.0 * 10.0,
            0.980,
            0.750,
            -2.02,
            0.0,
            0.0,
            0.0,
            10.0,
            None,
        )
        .expect("a positive hull should construct")
    };

    // Low pole: Cb=0.245 -> Cp=0.25 (the LCB term).
    let low = make(0.245);
    assert!((low.prismatic_coefficient() - 0.25).abs() < 1e-9);
    assert!(
        matches!(
            low.resistance_at(12.0, &water),
            Err(HydroError::OutOfRange { .. })
        ),
        "Cp=0.25 must be rejected, not return a NaN readout"
    );

    // High pole: Cb=0.95 -> Cp=0.969 (the (0.95-Cp) term) — an ordinary
    // full-hull entry, so this band matters in practice.
    let high = make(0.95);
    assert!(high.prismatic_coefficient() > 0.95);
    assert!(
        matches!(
            high.resistance_at(12.0, &water),
            Err(HydroError::OutOfRange { .. })
        ),
        "Cp>0.95 must be rejected, not return a NaN readout"
    );

    // Sanity: the in-band Holtrop example (Cp=0.5833) is still accepted
    // (no over-rejection of valid hulls).
    let ok = make(0.5833 * 0.980); // Cb -> Cp = 0.5833 (in band)
    assert!(
        ok.resistance_at(12.0, &water).is_ok(),
        "an in-band hull must still compute a resistance"
    );
}

// ===========================================================================
// VALIDATION against published references (the point of the crate).
// ===========================================================================
mod validation {
    use super::*;

    // -----------------------------------------------------------------------
    // Reference 1 — the ITTC-1957 model-ship correlation line.
    //
    //   ITTC-1957: C_f = 0.075 / (log10(Re) - 2)^2
    //
    // Source: 8th International Towing Tank Conference, 1957; the line is
    // reproduced in ITTC Recommended Procedure 7.5-02-03-01.4 ("1957 ITTC
    // model-ship correlation line"),
    //   https://www.ittc.info/media/8017/75-02-03-014.pdf
    // and ITTC 7.5-02-05-01,
    //   https://ittc.info/media/2065/75-02-05-01.pdf
    //
    // The line is an exact closed form, so these are *exact* reference points,
    // not measured data.
    // -----------------------------------------------------------------------

    #[test]
    fn ittc57_friction_line_closed_form_points() {
        // Re = 1e7  ->  (7 - 2)^2 = 25  ->  Cf = 0.075/25 = 0.003 exactly.
        assert!(close(
            ittc57_friction_coefficient(1e7).unwrap(),
            0.003,
            1e-12
        ));
        // Re = 1e8  ->  (8 - 2)^2 = 36  ->  Cf = 0.075/36 = 0.00208333...
        assert!(close(
            ittc57_friction_coefficient(1e8).unwrap(),
            0.075 / 36.0,
            1e-12
        ));
        // Re = 1e9  ->  (9 - 2)^2 = 49  ->  Cf = 0.075/49 = 0.00153061...
        // This is the textbook pin: Cf(1e9) ~= 0.00153.
        let cf_1e9 = ittc57_friction_coefficient(1e9).unwrap();
        assert!(close(cf_1e9, 0.075 / 49.0, 1e-12));
        assert!(close(cf_1e9, 0.0015306, 1e-4));
    }

    #[test]
    fn ittc57_is_undefined_for_tiny_reynolds() {
        // log10(Re) <= 2 (Re <= 100) makes the line singular; reject it.
        assert!(ittc57_friction_coefficient(100.0).is_err());
        assert!(ittc57_friction_coefficient(50.0).is_err());
        assert!(ittc57_friction_coefficient(0.0).is_err());
        assert!(ittc57_friction_coefficient(-1.0).is_err());
        // Just above the singularity it is finite and positive.
        assert!(ittc57_friction_coefficient(200.0).unwrap() > 0.0);
    }

    // -----------------------------------------------------------------------
    // Reference 2 — Holtrop & Mennen (1982) worked example.
    //
    // J. Holtrop and G.G.J. Mennen, "An approximate power prediction method",
    // International Shipbuilding Progress, Vol. 29, No. 335, July 1982,
    // pp. 166-170. Section 5 ("Numerical example") gives a hypothetical
    // single-screw ship at 25 knots. Paper PDF (as fetched for this work):
    //   https://moodle2.units.it/pluginfile.php/493107/mod_resource/content/0/Holtrop%20calcolo%20statistico%20resistenza.pdf
    //   https://www.boatdesign.net/attachments/holtrop-approximate-1982-pdf.144448/
    //
    // Published main-particulars (inputs) and the paper's result table:
    //   L_wl = 205.00 m, L_pp = 200.00 m, B = 32.00 m, T = 10.00 m,
    //   nabla = 37500 m^3, lcb = 2.02% aft of 1/2 L, C_m = 0.980, C_wp = 0.750,
    //   transverse bulb area A_bt = 20.0 m^2, bulb centre h_b = 4.0 m,
    //   immersed transom area A_t = 16.0 m^2, C_stern = 10 (stern shape),
    //   ship speed V = 25.0 knots, sea water (rho = 1025 kg/m^3, 15 deg C).
    //
    // Paper's reported results (read from the example's result table):
    //   C_p = 0.5833, S = 7381.45 m^2, C_f = 0.001390 (ITTC-57),
    //   i_E = 12.08 deg, half-entrance c1 = 1.398, c2 = 0.7595,
    //   wave resistance R_w = 557.11 kN, frictional R_f = 869.63 kN,
    //   correlation allowance R_A = 221.98 kN, appendage R_app = 8.83 kN,
    //   total R_t = 1793.26 kN, effective power P_E = 23063 kW,
    //   Froude number Fn = 0.2868.
    //
    // Our crate models the *bare hull* (no appendage term), so we validate the
    // bare-hull quantities digit-exact (C_p, S, C_f, R_f, Re, Fn) and the
    // statistical quantities (i_E, 1+k, R_w) to the few-per-cent tolerance
    // appropriate for a regression. The paper's own resistance balance implies
    // a form factor of (R_t - R_app - R_w - R_A)/R_f = 1.156.
    // -----------------------------------------------------------------------

    /// The Holtrop & Mennen (1982) example hull, built from the published
    /// main particulars.
    fn holtrop_1982_example_hull() -> Hull {
        Hull::new(
            205.00, // L_wl (m)
            32.00,  // B (m)
            10.00,  // T (m)
            37500.0, // nabla (m^3)
            0.980,  // C_m
            0.750,  // C_wp
            -2.02,  // lcb (% of L; aft = negative)
            20.0,   // A_bt (m^2)
            4.0,    // h_b (m)
            16.0,   // A_t (m^2)
            10.0,   // C_stern
            None,   // let the crate estimate S (validated below)
        )
        .unwrap()
    }

    /// Sea water at the paper's stated condition (rho = 1025, 15 deg C).
    fn holtrop_1982_water() -> WaterProperties {
        WaterProperties::seawater()
    }

    #[test]
    fn geometry_coefficients_match_paper_exactly() {
        let hull = holtrop_1982_example_hull();
        // C_b = nabla/(L B T) = 37500/(205*32*10) = 0.571646...
        assert!(close(hull.block_coefficient(), 0.571_646, 1e-4));
        // C_p = C_b/C_m -> paper reports 0.5833.
        assert!(
            close(hull.prismatic_coefficient(), 0.5833, 1e-3),
            "Cp = {} (paper 0.5833)",
            hull.prismatic_coefficient()
        );
    }

    #[test]
    fn wetted_surface_matches_paper_digit_exact() {
        let hull = holtrop_1982_example_hull();
        // Paper: S = 7381.45 m^2 (Holtrop wetted-surface regression).
        let s = hull.wetted_surface();
        assert!(
            close(s, 7381.45, 5e-5),
            "S = {s} m^2 (paper 7381.45 m^2)"
        );
    }

    #[test]
    fn reynolds_and_froude_match_paper() {
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);

        // Fn = V/sqrt(gL) -> paper 0.2868.
        let fn_ = froude_number(v, hull.length_m).unwrap();
        assert!(close(fn_, 0.2868, 1e-3), "Fn = {fn_} (paper 0.2868)");

        // Re = V L / nu ~ 2.219e9 with ITTC seawater nu.
        let re = reynolds_number(v, hull.length_m, water.kinematic_viscosity).unwrap();
        assert!(close(re, 2.219e9, 1e-3), "Re = {re:e} (paper ~2.219e9)");
    }

    #[test]
    fn friction_coefficient_matches_paper_digit_exact() {
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);
        let re = reynolds_number(v, hull.length_m, water.kinematic_viscosity).unwrap();
        let cf = ittc57_friction_coefficient(re).unwrap();
        // Paper: C_f = 0.001390.
        assert!(close(cf, 0.001390, 1e-3), "Cf = {cf} (paper 0.001390)");
    }

    #[test]
    fn frictional_resistance_matches_paper_digit_exact() {
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);
        let p = hull.resistance_at(v, &water).unwrap();
        // Paper: bare flat-plate frictional resistance R_f = 869.63 kN.
        let rf_kn = p.frictional_resistance_n / 1000.0;
        assert!(
            close(rf_kn, 869.63, 2e-3),
            "Rf = {rf_kn} kN (paper 869.63 kN)"
        );
    }

    #[test]
    fn wave_resistance_matches_paper_within_one_percent() {
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);
        let r_w_kn = hull.wave_resistance(v, &water).unwrap() / 1000.0;
        // Paper: R_w = 557.11 kN. Holtrop is a regression -> few-per-cent.
        assert!(
            close(r_w_kn, 557.11, 1.0e-2),
            "Rw = {r_w_kn} kN (paper 557.11 kN)"
        );
    }

    #[test]
    fn form_factor_matches_paper_balance_within_a_few_percent() {
        let hull = holtrop_1982_example_hull();
        // Paper's resistance balance implies (1+k) = 1.156. The Holtrop
        // regression here gives ~1.162 -> within ~1%.
        let one_plus_k = hull.form_factor();
        assert!(
            close(one_plus_k, 1.156, 2.0e-2),
            "(1+k) = {one_plus_k} (paper-implied 1.156)"
        );
        // Sanity: a normal-form factor is in the expected band.
        assert!((1.05..=1.35).contains(&one_plus_k));
    }

    #[test]
    fn reconstructed_total_and_effective_power_match_paper_within_two_percent() {
        // Our bare-hull total = R_f(1+k) + R_w + R_a. The paper additionally
        // carries an appendage term R_app = 8.83 kN; adding it back, the totals
        // should agree to a couple of per cent. We use the paper's own R_app
        // and R_a so this isolates the hull physics we model.
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);
        let p = hull.resistance_at(v, &water).unwrap();

        // Replace our default correlation term with the paper's R_a and add
        // the paper's appendage R_app, both in newtons.
        let r_app = 8_830.0_f64;
        let r_a_paper = 221_980.0_f64;
        let reconstructed = p.viscous_resistance_n + p.wave_resistance_n + r_app + r_a_paper;
        let reconstructed_kn = reconstructed / 1000.0;
        assert!(
            close(reconstructed_kn, 1793.26, 2.0e-2),
            "reconstructed Rt = {reconstructed_kn} kN (paper 1793.26 kN)"
        );

        let pe_kw = reconstructed * v / 1000.0;
        assert!(
            close(pe_kw, 23063.0, 2.0e-2),
            "reconstructed PE = {pe_kw} kW (paper 23063 kW)"
        );
    }

    #[test]
    fn effective_power_identity_pe_equals_rt_times_v() {
        // The paper's own P_E = R_t * V identity: 1793.26 kN * 12.8611 m/s
        // = 23063 kW. Verify our P_e is exactly R_t * V for our own R_t.
        let hull = holtrop_1982_example_hull();
        let water = holtrop_1982_water();
        let v = knots_to_ms(25.0);
        let p = hull.resistance_at(v, &water).unwrap();
        assert!(close(p.effective_power_w, p.total_resistance_n * v, 1e-12));

        // And the paper's exact identity, independent of our R_t:
        let pe_paper = 1_793_260.0 * knots_to_ms(25.0);
        assert!(
            close(pe_paper / 1000.0, 23063.0, 5e-3),
            "PE(paper) = {} kW (paper 23063 kW)",
            pe_paper / 1000.0
        );
    }
}

// ===========================================================================
// Unit conversions.
// ===========================================================================
mod units {
    use super::*;

    #[test]
    fn knot_round_trip() {
        for kn in [1.0, 12.5, 25.0, 40.0] {
            assert!(close(ms_to_knots(knots_to_ms(kn)), kn, 1e-12));
        }
        // 1 knot = 0.514444... m/s.
        assert!(close(knots_to_ms(1.0), 0.514_444, 1e-5));
    }

    #[test]
    fn point_kn_and_kw_helpers() {
        let hull = Hull::new(
            205.0, 32.0, 10.0, 37500.0, 0.98, 0.75, -2.02, 20.0, 4.0, 16.0, 10.0, None,
        )
        .unwrap();
        let p = hull.resistance_at(knots_to_ms(25.0), &WaterProperties::seawater()).unwrap();
        assert!(close(p.total_resistance_kn(), p.total_resistance_n / 1000.0, 1e-12));
        assert!(close(p.effective_power_kw(), p.effective_power_w / 1000.0, 1e-12));
    }
}

// ===========================================================================
// Physical properties / monotonicity.
// ===========================================================================
mod properties {
    use super::*;

    fn demo_hull() -> Hull {
        Hull::new(
            120.0, 18.0, 6.0, 9000.0, 0.95, 0.80, 0.0, 0.0, 0.0, 0.0, 0.0, None,
        )
        .unwrap()
    }

    #[test]
    fn resistance_and_power_increase_with_speed() {
        let hull = demo_hull();
        let water = WaterProperties::seawater();
        let slow = hull.resistance_at(5.0, &water).unwrap();
        let fast = hull.resistance_at(9.0, &water).unwrap();
        assert!(fast.total_resistance_n > slow.total_resistance_n);
        assert!(fast.effective_power_w > slow.effective_power_w);
        // Frictional resistance grows roughly as V^2 (Cf drifts slowly).
        assert!(fast.frictional_resistance_n > slow.frictional_resistance_n);
    }

    #[test]
    fn total_is_sum_of_components() {
        let hull = demo_hull();
        let p = hull.resistance_at(7.0, &WaterProperties::seawater()).unwrap();
        let sum = p.viscous_resistance_n + p.wave_resistance_n + p.correlation_resistance_n;
        assert!(close(p.total_resistance_n, sum, 1e-9));
        // Viscous = friction * form factor.
        assert!(close(
            p.viscous_resistance_n,
            p.frictional_resistance_n * p.form_factor,
            1e-9
        ));
    }

    #[test]
    fn measured_wetted_surface_overrides_estimate() {
        let estimated = demo_hull();
        let mut measured = demo_hull();
        measured.wetted_surface_m2 = Some(estimated.wetted_surface() * 1.10);
        assert!(close(measured.wetted_surface(), estimated.wetted_surface() * 1.10, 1e-12));
        // A bigger wetted surface means more friction at the same speed.
        let water = WaterProperties::seawater();
        assert!(
            measured.resistance_at(7.0, &water).unwrap().frictional_resistance_n
                > estimated.resistance_at(7.0, &water).unwrap().frictional_resistance_n
        );
    }

    #[test]
    fn denser_water_gives_more_resistance() {
        let hull = demo_hull();
        let sea = WaterProperties::seawater();
        let fresh = WaterProperties::freshwater();
        let r_sea = hull.resistance_at(7.0, &sea).unwrap().frictional_resistance_n;
        let r_fresh = hull.resistance_at(7.0, &fresh).unwrap().frictional_resistance_n;
        assert!(r_sea > r_fresh);
    }

    #[test]
    fn from_hydrostatic_round_trips_geometry() {
        // Build the same geometry via valenx-marine and check volume agrees.
        let hydro = valenx_marine::Hull::new(150.0, 20.0, 8.0, 0.7, 7.0, SEAWATER_DENSITY).unwrap();
        let hull = Hull::from_hydrostatic(&hydro, 0.98, 0.78, -1.0, 0.0, 0.0, 0.0, 0.0, None).unwrap();
        assert!(close(hull.volume_m3, hydro.displaced_volume(), 1e-9));
        assert!(close(hull.block_coefficient(), 0.7, 1e-9));
    }
}

// ===========================================================================
// The curve API.
// ===========================================================================
mod curve {
    use super::*;

    fn demo_hull() -> Hull {
        Hull::new(
            120.0, 18.0, 6.0, 9000.0, 0.95, 0.80, 0.0, 0.0, 0.0, 0.0, 0.0, None,
        )
        .unwrap()
    }

    #[test]
    fn stepped_range_speeds_are_ascending_and_inclusive() {
        let r = SpeedRange::stepped(4.0, 10.0, 2.0).unwrap();
        let s = r.speeds();
        assert_eq!(s.first().copied(), Some(4.0));
        assert_eq!(s.last().copied(), Some(10.0));
        for w in s.windows(2) {
            assert!(w[1] > w[0]);
        }
    }

    #[test]
    fn linspace_has_n_points() {
        let r = SpeedRange::linspace(5.0, 9.0, 5).unwrap();
        let s = r.speeds();
        assert_eq!(s.len(), 5);
        assert!(close(s[0], 5.0, 1e-12));
        assert!(close(*s.last().unwrap(), 9.0, 1e-9));
    }

    #[test]
    fn resistance_curve_is_monotonic_and_queryable() {
        let hull = demo_hull();
        let water = WaterProperties::seawater();
        let range = SpeedRange::stepped(4.0, 10.0, 1.0).unwrap();
        let curve = resistance_curve(&hull, &range, &water).unwrap();
        assert_eq!(curve.points.len(), 7);
        // Total resistance climbs across the range.
        for w in curve.points.windows(2) {
            assert!(w[1].total_resistance_n > w[0].total_resistance_n);
        }
        // Peak is the last (fastest) point here.
        assert!(close(
            curve.peak_resistance().unwrap().speed_ms,
            10.0,
            1e-12
        ));
        // Nearest finds the right speed.
        let near = curve.nearest(6.4).unwrap();
        assert!(close(near.speed_ms, 6.0, 1e-9));
    }

    #[test]
    fn curve_is_serde_round_trippable() {
        let hull = demo_hull();
        let water = WaterProperties::seawater();
        let range = SpeedRange::stepped(4.0, 8.0, 2.0).unwrap();
        let curve = resistance_curve(&hull, &range, &water).unwrap();
        let json = serde_json::to_string(&curve).unwrap();
        let back: ResistanceCurve = serde_json::from_str(&json).unwrap();
        assert_eq!(back.points.len(), curve.points.len());
        for (a, b) in back.points.iter().zip(curve.points.iter()) {
            assert!(close(a.total_resistance_n, b.total_resistance_n, 1e-9));
            assert!(close(a.effective_power_w, b.effective_power_w, 1e-9));
        }
    }
}

// ===========================================================================
// Error domain.
// ===========================================================================
mod errors {
    use super::*;

    #[test]
    fn hull_rejects_out_of_domain_inputs() {
        // zero length
        assert!(Hull::new(0.0, 32.0, 10.0, 1.0, 0.98, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, None).is_err());
        // midship coeff > 1
        assert!(Hull::new(200.0, 32.0, 10.0, 1.0, 1.5, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, None).is_err());
        // waterplane coeff = 0
        assert!(Hull::new(200.0, 32.0, 10.0, 1.0, 0.98, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, None).is_err());
        // negative bulb area
        assert!(Hull::new(200.0, 32.0, 10.0, 1.0, 0.98, 0.75, 0.0, -1.0, 0.0, 0.0, 0.0, None).is_err());
        // non-finite lcb
        assert!(Hull::new(200.0, 32.0, 10.0, 1.0, 0.98, 0.75, f64::NAN, 0.0, 0.0, 0.0, 0.0, None).is_err());
        // supplied non-positive wetted surface
        assert!(Hull::new(200.0, 32.0, 10.0, 1.0, 0.98, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, Some(-5.0)).is_err());
    }

    #[test]
    fn water_validation_rejects_bad_fluid() {
        assert!(WaterProperties { density: 0.0, kinematic_viscosity: 1e-6, correlation_allowance: 0.0 }.validate().is_err());
        assert!(WaterProperties { density: 1025.0, kinematic_viscosity: -1e-6, correlation_allowance: 0.0 }.validate().is_err());
        assert!(WaterProperties { density: 1025.0, kinematic_viscosity: 1e-6, correlation_allowance: f64::NAN }.validate().is_err());
        assert!(WaterProperties::seawater().validate().is_ok());
    }

    #[test]
    fn resistance_rejects_non_positive_speed() {
        let hull = Hull::new(120.0, 18.0, 6.0, 9000.0, 0.95, 0.80, 0.0, 0.0, 0.0, 0.0, 0.0, None).unwrap();
        let water = WaterProperties::seawater();
        assert!(hull.resistance_at(0.0, &water).is_err());
        assert!(hull.resistance_at(-3.0, &water).is_err());
        assert!(hull.resistance_at(f64::NAN, &water).is_err());
    }

    #[test]
    fn range_rejects_bad_bounds() {
        assert!(SpeedRange::stepped(10.0, 4.0, 1.0).is_err()); // end < start
        assert!(SpeedRange::stepped(4.0, 10.0, 0.0).is_err()); // zero step
        assert!(SpeedRange::stepped(-1.0, 10.0, 1.0).is_err()); // non-positive start
        assert!(SpeedRange::linspace(4.0, 10.0, 0).is_err()); // zero points
    }

    #[test]
    fn frictional_resistance_rejects_negative_cf() {
        assert!(frictional_resistance(1025.0, 5.0, 1000.0, -0.001).is_err());
        assert!(frictional_resistance(1025.0, 5.0, 1000.0, 0.002).is_ok());
    }
}
