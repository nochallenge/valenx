//! # valenx-osmosis — osmosis & fluid balance
//!
//! Closed-form textbook models for water movement across semipermeable
//! membranes: colligative osmotic pressure, transcapillary fluid
//! exchange, and tonicity classification. Pure arithmetic on `f64`; no
//! solvers, no external processes, no platform dependencies.
//!
//! ## What
//!
//! Three small, independent model families, each with validated
//! constructors and ground-truth-checked formulae.
//!
//! 1. Van't Hoff osmotic pressure ([`vant_hoff`]). The colligative law
//!    `Pi = i * c * R * T` for a dilute solution, the osmolarity
//!    `i * c`, results in atmospheres and pascals, and two inverses —
//!    the osmolarity from a measured pressure, and the solute **molar
//!    mass** (membrane osmometry) via
//!    [`vant_hoff::molar_mass_from_pressure_atm`]. See
//!    [`vant_hoff::Solution`].
//! 2. Starling filtration ([`starling`]). The transcapillary flux
//!    `Jv = Kf * ((Pc - Pi) - sigma * (pi_c - pi_i))`, decomposed into
//!    its net hydrostatic gradient, net oncotic gradient, net driving
//!    pressure, and a filtration / reabsorption / equilibrium direction
//!    classifier. See [`starling::StarlingParams`].
//! 3. Tonicity ([`tonicity`]). Classification of an external solution as
//!    hypotonic, isotonic, or hypertonic relative to a reference
//!    compartment by comparing effective (non-penetrating) osmolarities,
//!    plus the implied net water-movement direction. See
//!    [`tonicity::Tonicity`].
//!
//! ## Model
//!
//! The three families correspond to the standard physical-chemistry and
//! physiology relations.
//!
//! Osmotic pressure follows the van't Hoff ideal-dilute law, the
//! osmotic analogue of the ideal-gas law `Pi * V = n * R * T`. The unit
//! of `Pi` is set by the gas constant chosen: [`vant_hoff::R_L_ATM`]
//! with concentration in `mol/L` yields atmospheres,
//! [`vant_hoff::R_SI`] with concentration in `mol/m^3` yields pascals.
//! Pressure is therefore exactly proportional to both the particle
//! concentration `i * c` and the absolute temperature `T`. Because the
//! molar concentration is the mass concentration over the molar mass
//! (`c = rho / M`), a measured pressure also yields the solute's molar
//! mass, `M = i * rho * R * T / Pi` — the basis of membrane osmometry
//! for sizing macromolecules.
//!
//! Transcapillary flow follows Starling's equation: the difference
//! between the net hydrostatic pressure pushing fluid out of a
//! capillary and the reflection-coefficient-weighted net oncotic
//! pressure pulling it back in, scaled by the filtration coefficient.
//! Positive flux is filtration (out of the capillary); negative flux is
//! reabsorption.
//!
//! Tonicity compares the effective osmolarity (non-penetrating solutes
//! only) of an external solution against a cell-interior reference;
//! water moves toward the side of higher effective osmolarity, so equal
//! effective osmolarities produce no net movement.
//!
//! All pressures within a single computation share one unit (the crate
//! does not impose one; mmHg is conventional for Starling, atm or Pa for
//! van't Hoff), so callers must keep their inputs consistent.
//!
//! ## Honest scope
//!
//! Research / educational grade. Every formula here is a closed-form,
//! well-established textbook relation — the van't Hoff colligative law,
//! Starling's principle, and the osmolarity / tonicity comparison — and
//! the unit tests pin each against analytic and known physiological
//! values. It is deliberately **not** a clinical, medical, or
//! production engineering tool, and it makes the standard simplifying
//! assumptions:
//!
//! The van't Hoff law is the **ideal-dilute** limit. It assumes an
//! osmotic coefficient of exactly `1`; real electrolytes deviate (the
//! effective van't Hoff factor of NaCl in plasma is closer to `1.85`
//! than the ideal `2`, and concentrated solutions need an activity
//! correction). Pass an empirical `i` if you want to absorb that.
//!
//! The Starling model is the **steady, single-segment** form. It uses
//! constant lumped parameters and omits the dynamic revision of
//! Starling's principle (the sub-glycocalyx oncotic pressure / "no
//! reabsorption in steady state" picture), lymphatic return, and any
//! spatial variation of `Pc` along the capillary.
//!
//! Tonicity is computed from **caller-supplied effective osmolarities**.
//! The crate does not itself decide which solutes penetrate the
//! membrane; distinguishing tonicity from osmolarity (the urea /
//! cryoprotectant case) is the caller's responsibility, as documented
//! on [`tonicity`].
//!
//! ## Errors
//!
//! Every fallible constructor and computation returns
//! [`Result<_, OsmosisError>`](error::OsmosisError). The error exposes
//! stable [`code`](error::OsmosisError::code) and
//! [`category`](error::OsmosisError::category) accessors for telemetry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod starling;
pub mod tonicity;
pub mod vant_hoff;

pub use error::{ErrorCategory, OsmosisError};
pub use starling::{FluxDirection, StarlingParams};
pub use tonicity::{classify as classify_tonicity, Tonicity, WaterMovement};
pub use vant_hoff::{
    molar_mass_from_pressure_atm, osmolarity_from_pressure_atm, osmotic_pressure_atm, Solution,
    ATM_IN_PA, CELSIUS_ZERO_K, R_L_ATM, R_SI,
};

#[cfg(test)]
mod integration_tests {
    //! Cross-module checks tying the three families together on a single
    //! physiological scenario.

    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn isotonic_saline_has_matching_osmotic_pressure_and_no_net_water() {
        // Normal saline: ~0.300 osmol/L effective osmolarity at 37 C.
        // Model it as 0.150 mol/L NaCl with i = 2 -> 0.300 osmol/L.
        let saline = Solution::new(0.150, 2.0, 310.15).unwrap();
        let cell_osmolarity = 0.300;

        // Osmolarity matches the cell reference.
        assert!((saline.osmolarity_osmol_per_l() - cell_osmolarity).abs() < EPS);

        // Therefore the tonicity is isotonic and water does not move.
        let t = classify_tonicity(saline.osmolarity_osmol_per_l(), cell_osmolarity).unwrap();
        assert_eq!(t, Tonicity::Isotonic);
        assert_eq!(t.water_movement(), WaterMovement::NoNet);

        // And its colligative osmotic pressure lands in the textbook
        // plasma band (~7.6-7.7 atm at body temperature).
        let pi = saline.osmotic_pressure_atm();
        assert!((7.0..8.5).contains(&pi), "pi = {pi}");
    }

    #[test]
    fn capillary_filtration_balances_over_a_cycle() {
        // Arteriolar end filters (+14 mmHg) and venular end reabsorbs
        // (-6 mmHg) with the same membrane; the net is outward filtration
        // (the small surplus the lymphatics return in vivo).
        let arterial = StarlingParams::new(35.0, -2.0, 28.0, 5.0, 1.0, 1.0).unwrap();
        let venular = StarlingParams::new(15.0, -2.0, 28.0, 5.0, 1.0, 1.0).unwrap();

        assert_eq!(arterial.direction(), FluxDirection::Filtration);
        assert_eq!(venular.direction(), FluxDirection::Reabsorption);

        let net = arterial.net_filtration() + venular.net_filtration();
        // +14 + (-6) = +8 mmHg worth of net filtration (Kf = 1).
        assert!((net - 8.0).abs() < EPS, "net = {net}");
    }
}
