//! # valenx-combustion
//!
//! Combustion stoichiometry for pure hydrocarbon fuels `CxHy` burned in
//! air.
//!
//! ## What
//!
//! Closed-form thermochemistry of complete hydrocarbon combustion:
//!
//! - [`fuel::Fuel`] — a validated `CxHy` fuel with molar mass and the
//!   stoichiometric O2 demand `a = x + y/4`. Named constructors cover
//!   [`Fuel::methane`], [`Fuel::propane`], [`Fuel::octane`].
//! - [`stoich::afr_stoich_mass`] / [`stoich::afr_stoich_molar`] — the
//!   stoichiometric air-fuel ratio on a mass and a molar basis.
//! - [`stoich::equivalence_ratio`] / [`stoich::afr_from_phi`] —
//!   `phi = AFR_stoich / AFR_actual` and its inverse.
//! - [`stoich::percent_theoretical_air`] /
//!   [`stoich::percent_excess_air`] — the engineering air-side
//!   descriptors `100/phi` and `100(1 - phi)/phi`.
//! - [`stoich::classify`] — lean (`phi < 1`) / stoichiometric
//!   (`phi == 1`) / rich (`phi > 1`).
//! - [`stoich::product_moles`] — the balanced product mole table
//!   (CO2, H2O, N2, excess O2) per mole of fuel.
//! - [`flame::adiabatic_flame_temperature`] — a constant-pressure
//!   mean-cp energy-balance flame-temperature estimate.
//!
//! ## Model
//!
//! For one mole of fuel with `a = x + y/4` moles of O2 (carrying
//! `3.76 a` moles of inert N2 in the textbook dry-air model
//! `O2 + 3.76 N2`):
//!
//! ```text
//! CxHy + a (O2 + 3.76 N2) -> x CO2 + (y/2) H2O + 3.76 a N2
//! ```
//!
//! The stoichiometric mass air-fuel ratio is
//! `a (M_O2 + 3.76 M_N2) / M_fuel`; the equivalence ratio rescales the
//! actual air supply. The same mixture is equivalently described by the
//! percent theoretical air `100/phi` or percent excess air
//! `100(1 - phi)/phi` (negative — an air deficiency — when rich). The adiabatic flame temperature comes from a
//! first-law balance, `T_ad = T_in + LHV * M_fuel / (n_products *
//! cp_molar)`, where all of the released chemical energy heats the
//! product gas through a single constant mean molar heat capacity.
//!
//! ## Honest scope
//!
//! Research and educational grade. These are textbook closed-form
//! complete-combustion relations plus a single-mean-cp energy balance —
//! NOT a clinical/medical/production engineering tool. The flame-
//! temperature model ignores chemical dissociation, equilibrium
//! product speciation (CO, H2, NO, OH, soot), the strong temperature
//! dependence of cp, and finite-rate kinetics, so its temperatures are
//! first-order estimates that run several hundred kelvin high relative
//! to a real flame. The product balance is only defined for the lean
//! and stoichiometric regimes (`phi <= 1`); the rich regime would
//! require equilibrium chemistry and is reported as unsupported.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod flame;
pub mod fuel;
pub mod stoich;

pub use error::{CombustionError, ErrorCategory};
pub use flame::{
    adiabatic_flame_temperature, heat_release_per_mole_fuel, CP_MOLAR_DEFAULT, LHV_METHANE,
    LHV_OCTANE, LHV_PROPANE, T_REF_K,
};
pub use fuel::Fuel;
pub use stoich::{
    afr_from_phi, afr_stoich_mass, afr_stoich_molar, air_mean_molar_mass, classify,
    equivalence_ratio, percent_excess_air, percent_theoretical_air, product_moles, Mixture,
    ProductMoles,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for tight numeric checks.
    const EPS: f64 = 1e-9;

    // --- Fuel construction and validation -------------------------------

    #[test]
    fn fuel_rejects_zero_carbon() {
        let err = Fuel::new(0, 4).unwrap_err();
        assert_eq!(err.code(), "combustion.invalid_fuel");
        assert_eq!(err.category(), ErrorCategory::Input);
    }

    #[test]
    fn fuel_rejects_zero_hydrogen() {
        let err = Fuel::new(1, 0).unwrap_err();
        assert_eq!(err.code(), "combustion.invalid_fuel");
    }

    #[test]
    fn fuel_accepts_valid_formula() {
        let f = Fuel::new(3, 8).unwrap();
        assert_eq!(f.carbon, 3);
        assert_eq!(f.hydrogen, 8);
        assert_eq!(f, Fuel::propane());
    }

    #[test]
    fn methane_molar_mass_is_16() {
        // 1*12.011 + 4*1.008 = 16.043 g/mol.
        let m = Fuel::methane().molar_mass();
        assert!((m - 16.043).abs() < 1e-6, "methane molar mass {m}");
    }

    #[test]
    fn octane_molar_mass_is_114() {
        // 8*12.011 + 18*1.008 = 114.232 g/mol.
        let m = Fuel::octane().molar_mass();
        assert!((m - 114.232).abs() < 1e-6, "octane molar mass {m}");
    }

    #[test]
    fn stoich_o2_demand_methane() {
        // a = 1 + 4/4 = 2.
        let a = Fuel::methane().stoich_o2_moles();
        assert!((a - 2.0).abs() < EPS, "O2 demand {a}");
    }

    #[test]
    fn stoich_o2_demand_octane() {
        // a = 8 + 18/4 = 12.5.
        let a = Fuel::octane().stoich_o2_moles();
        assert!((a - 12.5).abs() < EPS, "O2 demand {a}");
    }

    // --- Stoichiometric AFR (ground truth) ------------------------------

    #[test]
    fn methane_stoich_afr_mass_is_about_17() {
        // Textbook value for CH4 is ~17.1-17.2 kg air / kg fuel.
        let afr = afr_stoich_mass(&Fuel::methane());
        assert!((afr - 17.12).abs() < 0.05, "CH4 AFR {afr}");
    }

    #[test]
    fn octane_stoich_afr_mass_is_about_15() {
        // Gasoline / iso-octane textbook value ~15.0-15.1.
        let afr = afr_stoich_mass(&Fuel::octane());
        assert!((afr - 15.03).abs() < 0.1, "C8H18 AFR {afr}");
    }

    #[test]
    fn propane_stoich_afr_mass_is_about_15_6() {
        // Propane textbook value ~15.6-15.7.
        let afr = afr_stoich_mass(&Fuel::propane());
        assert!((afr - 15.6).abs() < 0.15, "C3H8 AFR {afr}");
    }

    #[test]
    fn methane_stoich_afr_molar() {
        // a*(1+3.76) = 2 * 4.76 = 9.52 mol air / mol fuel.
        let afr = afr_stoich_molar(&Fuel::methane());
        assert!((afr - 9.52).abs() < EPS, "CH4 molar AFR {afr}");
    }

    // --- Equivalence ratio ----------------------------------------------

    #[test]
    fn phi_is_one_at_stoichiometric() {
        let fuel = Fuel::methane();
        let afr_s = afr_stoich_mass(&fuel);
        let phi = equivalence_ratio(&fuel, afr_s).unwrap();
        assert!((phi - 1.0).abs() < EPS, "stoich phi {phi}");
    }

    #[test]
    fn phi_greater_than_one_is_rich() {
        // Less air than stoichiometric => rich => phi > 1.
        let fuel = Fuel::methane();
        let afr_s = afr_stoich_mass(&fuel);
        let phi = equivalence_ratio(&fuel, afr_s * 0.5).unwrap();
        assert!(phi > 1.0, "rich phi {phi}");
        assert!((phi - 2.0).abs() < EPS, "rich phi {phi}");
    }

    #[test]
    fn phi_less_than_one_is_lean() {
        // More air than stoichiometric => lean => phi < 1.
        let fuel = Fuel::methane();
        let afr_s = afr_stoich_mass(&fuel);
        let phi = equivalence_ratio(&fuel, afr_s * 2.0).unwrap();
        assert!(phi < 1.0, "lean phi {phi}");
        assert!((phi - 0.5).abs() < EPS, "lean phi {phi}");
    }

    #[test]
    fn afr_from_phi_round_trips() {
        let fuel = Fuel::octane();
        let phi_target = 0.8;
        let afr = afr_from_phi(&fuel, phi_target).unwrap();
        let phi_back = equivalence_ratio(&fuel, afr).unwrap();
        assert!(
            (phi_back - phi_target).abs() < EPS,
            "round-trip phi {phi_back}"
        );
    }

    #[test]
    fn equivalence_ratio_rejects_nonpositive_afr() {
        let err = equivalence_ratio(&Fuel::methane(), 0.0).unwrap_err();
        assert_eq!(err.code(), "combustion.bad_parameter");
        assert_eq!(err.category(), ErrorCategory::Config);
    }

    #[test]
    fn equivalence_ratio_rejects_nan_and_inf() {
        // The positive-finite guard must catch NaN and infinity, not just
        // negatives and zero.
        assert!(equivalence_ratio(&Fuel::methane(), f64::NAN).is_err());
        assert!(equivalence_ratio(&Fuel::methane(), f64::INFINITY).is_err());
        assert!(equivalence_ratio(&Fuel::methane(), -5.0).is_err());
    }

    // --- Percent theoretical / excess air -------------------------------

    #[test]
    fn theoretical_air_is_hundred_at_stoichiometric() {
        assert!((percent_theoretical_air(1.0).unwrap() - 100.0).abs() < EPS);
        assert!(percent_excess_air(1.0).unwrap().abs() < EPS);
    }

    #[test]
    fn theoretical_and_excess_air_closed_form() {
        // phi = 0.5 (double the air): 200% theoretical, 100% excess.
        assert!((percent_theoretical_air(0.5).unwrap() - 200.0).abs() < EPS);
        assert!((percent_excess_air(0.5).unwrap() - 100.0).abs() < EPS);
        // phi = 0.8: 125% theoretical, 25% excess.
        assert!((percent_theoretical_air(0.8).unwrap() - 125.0).abs() < 1e-9);
        assert!((percent_excess_air(0.8).unwrap() - 25.0).abs() < 1e-9);
    }

    #[test]
    fn theoretical_air_is_excess_plus_hundred() {
        // Exact identity at any phi (lean, stoichiometric, and rich).
        for &phi in &[0.4, 0.7, 1.0, 1.5, 2.0] {
            let th = percent_theoretical_air(phi).unwrap();
            let ex = percent_excess_air(phi).unwrap();
            assert!((th - (ex + 100.0)).abs() < 1e-9, "identity at phi={phi}");
        }
    }

    #[test]
    fn rich_mixture_has_air_deficiency() {
        // phi = 2 (half the air): 50% theoretical, -50% excess.
        assert!((percent_theoretical_air(2.0).unwrap() - 50.0).abs() < EPS);
        assert!((percent_excess_air(2.0).unwrap() + 50.0).abs() < EPS);
    }

    #[test]
    fn excess_air_consistent_with_afr() {
        // AFR_actual = AFR_stoich * (1 + excess/100): tie the excess-air
        // representation back to the actual air-fuel ratio.
        let fuel = Fuel::octane();
        let phi = 0.75;
        let ex = percent_excess_air(phi).unwrap();
        let afr_actual = afr_from_phi(&fuel, phi).unwrap();
        let afr_stoich = afr_stoich_mass(&fuel);
        assert!(
            (afr_actual - afr_stoich * (1.0 + ex / 100.0)).abs() < 1e-9,
            "AFR vs excess-air mismatch"
        );
    }

    #[test]
    fn more_excess_air_means_leaner() {
        // Lower phi (leaner) -> more excess air.
        let lean = percent_excess_air(0.5).unwrap();
        let near = percent_excess_air(0.9).unwrap();
        assert!(
            lean > near,
            "leaner should have more excess air: {lean} vs {near}"
        );
    }

    #[test]
    fn air_descriptors_reject_nonpositive_phi() {
        assert!(percent_theoretical_air(0.0).is_err());
        assert!(percent_excess_air(-1.0).is_err());
        assert!(percent_theoretical_air(f64::NAN).is_err());
        assert!(percent_excess_air(f64::INFINITY).is_err());
    }

    // --- Classification -------------------------------------------------

    #[test]
    fn classify_buckets() {
        assert_eq!(classify(0.7, 1e-6).unwrap(), Mixture::Lean);
        assert_eq!(classify(1.0, 1e-6).unwrap(), Mixture::Stoichiometric);
        assert_eq!(classify(1.3, 1e-6).unwrap(), Mixture::Rich);
    }

    #[test]
    fn classify_tolerance_band() {
        // Within tolerance of 1.0 reads as stoichiometric.
        assert_eq!(classify(1.02, 0.05).unwrap(), Mixture::Stoichiometric);
        assert_eq!(classify(0.98, 0.05).unwrap(), Mixture::Stoichiometric);
        // Outside it tips lean / rich.
        assert_eq!(classify(1.10, 0.05).unwrap(), Mixture::Rich);
    }

    #[test]
    fn classify_rejects_bad_inputs() {
        assert!(classify(0.0, 1e-6).is_err());
        assert!(classify(-1.0, 1e-6).is_err());
        assert!(classify(1.0, -1e-6).is_err());
    }

    // --- Product mole balance -------------------------------------------

    #[test]
    fn product_balance_methane_stoich() {
        // CH4 + 2(O2 + 3.76 N2) -> CO2 + 2 H2O + 7.52 N2, no excess O2.
        let p = product_moles(&Fuel::methane(), 1.0).unwrap();
        assert!((p.co2 - 1.0).abs() < EPS, "co2 {}", p.co2);
        assert!((p.h2o - 2.0).abs() < EPS, "h2o {}", p.h2o);
        assert!((p.n2 - 7.52).abs() < EPS, "n2 {}", p.n2);
        assert!(p.excess_o2.abs() < EPS, "excess o2 {}", p.excess_o2);
        assert!((p.total() - 10.52).abs() < EPS, "total {}", p.total());
    }

    #[test]
    fn product_balance_conserves_carbon_and_hydrogen() {
        // Carbon: co2 == x. Hydrogen: 2*h2o == y. Holds at any phi<=1.
        let fuel = Fuel::octane();
        for &phi in &[0.5, 0.75, 1.0] {
            let p = product_moles(&fuel, phi).unwrap();
            assert!(
                (p.co2 - fuel.carbon as f64).abs() < EPS,
                "carbon at phi={phi}"
            );
            assert!(
                (2.0 * p.h2o - fuel.hydrogen as f64).abs() < EPS,
                "hydrogen at phi={phi}"
            );
        }
    }

    #[test]
    fn product_balance_conserves_oxygen_atoms() {
        // O atoms supplied = 2 * O2_supplied. O atoms in products =
        // 2*co2 + h2o + 2*excess_o2. Must match exactly.
        let fuel = Fuel::propane();
        let phi = 0.6;
        let p = product_moles(&fuel, phi).unwrap();
        let a_supplied = fuel.stoich_o2_moles() / phi;
        let o_in = 2.0 * a_supplied;
        let o_out = 2.0 * p.co2 + p.h2o + 2.0 * p.excess_o2;
        assert!(
            (o_in - o_out).abs() < 1e-9,
            "O balance in={o_in} out={o_out}"
        );
    }

    #[test]
    fn lean_mixture_has_excess_oxygen() {
        let p = product_moles(&Fuel::methane(), 0.5).unwrap();
        // phi=0.5 => double the air => excess O2 = a_stoich*(1/phi - 1)
        // = 2*(2-1) = 2.
        assert!((p.excess_o2 - 2.0).abs() < EPS, "excess o2 {}", p.excess_o2);
    }

    #[test]
    fn rich_product_balance_is_unsupported() {
        let err = product_moles(&Fuel::methane(), 1.5).unwrap_err();
        assert_eq!(err.code(), "combustion.bad_parameter");
    }

    // --- Adiabatic flame temperature ------------------------------------

    #[test]
    fn heat_release_methane_is_about_802_kj_per_mol() {
        // LHV 50 MJ/kg * 0.016043 kg/mol = 802.15 kJ/mol — matches the
        // textbook methane LHV of ~802 kJ/mol.
        let q = heat_release_per_mole_fuel(&Fuel::methane(), LHV_METHANE).unwrap();
        assert!((q - 802_150.0).abs() < 1.0, "Q {q}");
    }

    #[test]
    fn methane_flame_temperature_is_physical() {
        // No-dissociation mean-cp estimate lands ~2000-2500 K for the
        // stoichiometric methane/air flame.
        let t = adiabatic_flame_temperature(
            &Fuel::methane(),
            1.0,
            LHV_METHANE,
            CP_MOLAR_DEFAULT,
            T_REF_K,
        )
        .unwrap();
        assert!((2000.0..2600.0).contains(&t), "flame T {t}");
    }

    #[test]
    fn flame_temperature_peaks_near_stoichiometric() {
        // Adding excess air (leaner) dilutes the heat over more product
        // moles, so the lean flame is cooler than the stoichiometric one.
        let fuel = Fuel::methane();
        let t_stoich =
            adiabatic_flame_temperature(&fuel, 1.0, LHV_METHANE, CP_MOLAR_DEFAULT, T_REF_K)
                .unwrap();
        let t_lean =
            adiabatic_flame_temperature(&fuel, 0.6, LHV_METHANE, CP_MOLAR_DEFAULT, T_REF_K)
                .unwrap();
        assert!(
            t_stoich > t_lean,
            "stoich {t_stoich} should exceed lean {t_lean}"
        );
    }

    #[test]
    fn flame_temperature_rises_with_heating_value() {
        // Same fuel geometry, higher LHV => more heat => hotter flame.
        let fuel = Fuel::methane();
        let t_lo =
            adiabatic_flame_temperature(&fuel, 1.0, 40.0e6, CP_MOLAR_DEFAULT, T_REF_K).unwrap();
        let t_hi =
            adiabatic_flame_temperature(&fuel, 1.0, 55.0e6, CP_MOLAR_DEFAULT, T_REF_K).unwrap();
        assert!(t_hi > t_lo, "higher LHV hotter: {t_hi} vs {t_lo}");
    }

    #[test]
    fn flame_temperature_balance_is_exact() {
        // Reconstruct T_ad by hand from the closed-form balance.
        let fuel = Fuel::propane();
        let phi = 0.9;
        let cp = 42.0;
        let t_in = 300.0;
        let q = heat_release_per_mole_fuel(&fuel, LHV_PROPANE).unwrap();
        let n = product_moles(&fuel, phi).unwrap().total();
        let expected = t_in + q / (n * cp);
        let got = adiabatic_flame_temperature(&fuel, phi, LHV_PROPANE, cp, t_in).unwrap();
        assert!(
            (got - expected).abs() < 1e-6,
            "got {got} expected {expected}"
        );
    }

    #[test]
    fn flame_temperature_rejects_bad_cp_and_tin() {
        let fuel = Fuel::methane();
        assert!(adiabatic_flame_temperature(&fuel, 1.0, LHV_METHANE, 0.0, T_REF_K).is_err());
        assert!(
            adiabatic_flame_temperature(&fuel, 1.0, LHV_METHANE, CP_MOLAR_DEFAULT, 0.0).is_err()
        );
    }

    #[test]
    fn heat_release_rejects_nonpositive_lhv() {
        let err = heat_release_per_mole_fuel(&Fuel::methane(), -1.0).unwrap_err();
        assert_eq!(err.code(), "combustion.bad_parameter");
    }

    // --- Air constant ----------------------------------------------------

    #[test]
    fn air_mean_molar_mass_is_about_28_8() {
        // (31.998 + 3.76*28.013)/4.76 ~= 28.85 g/mol — close to the
        // accepted 28.97 for real dry air (the 3.76 model is slightly
        // light because it lumps argon into N2).
        let m = air_mean_molar_mass();
        assert!((m - 28.85).abs() < 0.1, "air molar mass {m}");
    }
}
