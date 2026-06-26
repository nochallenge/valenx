//! Validation of `valenx-thermo` against named ground-truth values.
//!
//! Sources:
//! * Ideal-gas limit `Z → 1` as `P → 0` — thermodynamic identity.
//! * CO₂ compressibility — NIST Chemistry WebBook (Span–Wagner reference EoS).
//! * Saturation pressures — NIST WebBook / DIPPR experimental values.
//!
//! Cubic EoS are engineering correlations, so the tolerances reflect their
//! known accuracy (a few percent on `Psat`, sub-percent on near-ideal `Z`).

use valenx_thermo::{saturation_pressure, Eos, EosModel, Fluid};

/// Relative error helper.
fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

#[test]
fn ideal_gas_limit_all_fluids() {
    // For every library fluid, the vapor-root Z must approach 1 as P → 0.
    for fluid in [
        Fluid::carbon_dioxide(),
        Fluid::nitrogen(),
        Fluid::methane(),
        Fluid::propane(),
    ] {
        for model in [EosModel::PengRobinson, EosModel::Srk] {
            let eos = Eos::new(fluid, model);
            let t = 1.5 * fluid.tc; // supercritical, single gas phase
            let z = eos.z_roots(t, 1.0).unwrap().vapor; // 1 Pa ≈ ideal
            assert!(
                rel(z, 1.0) < 1e-4,
                "{} {:?}: Z(P→0) = {z}",
                fluid.name,
                model
            );
        }
    }
}

#[test]
fn co2_compressibility_low_pressure() {
    // CO₂ gas at 313.15 K, 0.2 MPa is in the regime where a cubic EoS is
    // accurate. From the experimental second virial coefficient
    // B(313 K) ≈ −123 cm³/mol (Dymond & Smith, *The Virial Coefficients of
    // Pure Gases and Mixtures*), the leading virial estimate is
    //   Z ≈ 1 + B·P/(RT) ≈ 1 + (−123e-6)(0.2e6)/(8.314·313.15) ≈ 0.9905.
    // PR must reproduce this slightly-below-unity real-gas Z.
    let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
    let z = eos.z_roots(313.15, 0.2e6).unwrap().vapor;
    assert!(
        rel(z, 0.9905) < 0.002,
        "PR CO2 Z(313.15 K, 0.2 MPa) = {z}, virial-from-B ≈ 0.9905"
    );
}

#[test]
fn co2_second_virial_coefficient_vs_experiment() {
    // The low-density limit of a cubic EoS fixes its second virial coefficient
    //   B(T) = lim_{P→0} (Z − 1) R T / P.
    // For CO₂ at 313.15 K the experimental value is B ≈ −123 cm³/mol
    // (Dymond & Smith). Peng–Robinson should land within a few cm³/mol.
    use valenx_thermo::GAS_CONSTANT;
    let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
    let t = 313.15;
    let p = 1.0e5; // 0.1 MPa — low enough that the linear-in-P term dominates
    let z = eos.z_roots(t, p).unwrap().vapor;
    let b_cm3 = (z - 1.0) * GAS_CONSTANT * t / p * 1.0e6; // m³/mol → cm³/mol
    assert!(
        (b_cm3 - (-123.0)).abs() < 5.0,
        "PR CO2 B(313.15 K) = {b_cm3:.1} cm³/mol, exp ≈ -123"
    );
}

#[test]
fn nitrogen_compressibility_near_ideal() {
    // N₂ at 300 K, 1 MPa is nearly ideal: NIST Z ≈ 0.9989. PR within 0.5 %.
    let eos = Eos::new(Fluid::nitrogen(), EosModel::PengRobinson);
    let z = eos.z_roots(300.0, 1.0e6).unwrap().vapor;
    assert!(
        rel(z, 0.9989) < 0.005,
        "PR N2 Z(300 K, 1 MPa) = {z}, NIST ≈ 0.9989"
    );
}

#[test]
fn co2_saturation_pressure_vs_experiment() {
    // CO₂ at 273.15 K: experimental Psat = 3.4851 MPa (NIST).
    let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
    let psat = saturation_pressure(&eos, 273.15).unwrap();
    assert!(
        rel(psat, 3.4851e6) < 0.03,
        "PR CO2 Psat(273.15 K) = {psat:.4e}, exp ≈ 3.4851e6"
    );
}

#[test]
fn nitrogen_saturation_pressure_vs_experiment() {
    // N₂ at 100 K: experimental Psat ≈ 0.7773 MPa (NIST).
    let eos = Eos::new(Fluid::nitrogen(), EosModel::PengRobinson);
    let psat = saturation_pressure(&eos, 100.0).unwrap();
    assert!(
        rel(psat, 0.7773e6) < 0.05,
        "PR N2 Psat(100 K) = {psat:.4e}, exp ≈ 0.7773e6"
    );
}

#[test]
fn propane_saturation_pressure_vs_experiment() {
    // Propane at 300 K: experimental Psat ≈ 0.9979 MPa (NIST).
    // Check both PR and SRK land within a few percent.
    for model in [EosModel::PengRobinson, EosModel::Srk] {
        let eos = Eos::new(Fluid::propane(), model);
        let psat = saturation_pressure(&eos, 300.0).unwrap();
        assert!(
            rel(psat, 0.9979e6) < 0.04,
            "{model:?} C3H8 Psat(300 K) = {psat:.4e}, exp ≈ 0.9979e6"
        );
    }
}

#[test]
fn liquid_root_is_denser_than_vapor() {
    // Inside the two-phase window the liquid Z must be well below the vapor Z.
    let eos = Eos::new(Fluid::carbon_dioxide(), EosModel::PengRobinson);
    let psat = saturation_pressure(&eos, 273.15).unwrap();
    let roots = eos.z_roots(273.15, psat).unwrap();
    assert!(roots.three_real, "expected three real roots at saturation");
    assert!(
        roots.liquid < roots.vapor,
        "liquid Z {} should be < vapor Z {}",
        roots.liquid,
        roots.vapor
    );
}

#[test]
fn bad_input_fails_loud() {
    let eos = Eos::new(Fluid::methane(), EosModel::PengRobinson);
    assert!(eos.z_roots(-1.0, 1e5).is_err()); // negative T
    assert!(eos.z_roots(300.0, -1e5).is_err()); // negative P
    assert!(eos.pressure(300.0, -1e-3).is_err()); // negative V
    assert!(saturation_pressure(&eos, -5.0).is_err()); // negative T
    assert!(Fluid::new("x", -1.0, 1e6, 0.1).is_err()); // bad Tc
    assert!(Fluid::new("x", 300.0, -1.0, 0.1).is_err()); // bad Pc
}
