//! # valenx-heat-transfer — steady 1D heat transfer
//!
//! A small, dependency-light toolkit for the classic one-dimensional,
//! steady-state heat-transfer calculations taught in every introductory
//! heat-transfer course: conduction and convection thermal resistances,
//! the series/parallel thermal-resistance circuit, the straight-fin
//! heat rate and efficiency, and the linear conduction temperature
//! profile through a slab.
//!
//! ## What
//!
//! - **Conduction** ([`conduction`]) — the plane-wall conductive
//!   resistance `R = L / (k·A)`, the steady heat rate `Q = ΔT / R`, and
//!   the linear temperature profile `T(x)` through the wall via
//!   [`PlaneWall`].
//! - **Convection** ([`convection`]) — the surface-convection
//!   resistance `R = 1 / (h·A)` and Newton's-law-of-cooling heat rate
//!   via [`ConvectiveSurface`].
//! - **Resistance networks** ([`network`]) — series (`R = ΣRᵢ`) and
//!   parallel (`1/R = Σ1/Rᵢ`) combination of thermal resistors as a
//!   reducible [`ResistanceNetwork`] tree, plus the circuit heat rate
//!   `Q = ΔT / R_total`.
//! - **Fins** ([`fin`]) — the straight-fin base heat rate
//!   `q = √(h·P·k·A_c)·θ_b·tanh(m·L)` and the fin efficiency
//!   `η = tanh(mL)/(mL) ∈ (0, 1]` for adiabatic- and convective-tip
//!   boundary conditions via [`Fin`].
//!
//! ## Model
//!
//! Everything here is the **thermal–electrical analogue**: heat rate
//! `Q` plays the role of current and temperature difference `ΔT` the
//! role of voltage, so a steady 1D heat path reduces to
//!
//! ```text
//! Q = ΔT / R_total
//! ```
//!
//! Conductive layers in the direction of heat flow add their
//! resistances (series); alternative parallel paths add their
//! conductances. Fins are handled with the standard 1D fin equation
//! `d²θ/dx² = m²θ`, `m = √(hP / kA_c)`, whose closed-form solution
//! gives the base heat rate and the `tanh(mL)/(mL)` efficiency. All
//! formulae assume constant properties, a uniform film coefficient, and
//! genuine one-dimensionality.
//!
//! ## Honest scope
//!
//! Research/educational grade. These are well-established **textbook
//! closed-form / lumped-resistance** models (as in Incropera,
//! *Fundamentals of Heat and Mass Transfer*, or Cengel, *Heat and Mass
//! Transfer*) for steady, 1D, constant-property situations. They do
//! **not** model transient response, multidimensional conduction,
//! radiation, temperature-dependent properties, contact resistance, or
//! convection-coefficient correlations — and the crate is emphatically
//! **NOT a clinical/medical or production thermal-engineering tool**.
//! Use it to learn, prototype and sanity-check, not to certify
//! hardware. Each module documents its own assumptions inline.
//!
//! Every fallible call returns [`Result<_, HeatTransferError>`]; the
//! error type exposes stable [`code`](HeatTransferError::code) and
//! [`category`](HeatTransferError::category) accessors for telemetry.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod conduction;
pub mod convection;
pub mod error;
pub mod fin;
pub mod network;

pub use conduction::{conduction_resistance, PlaneWall};
pub use convection::{convection_resistance, ConvectiveSurface};
pub use error::{ErrorCategory, HeatTransferError, Result};
pub use fin::{Fin, TipCondition};
pub use network::{parallel_resistance, series_resistance, ResistanceNetwork};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// End-to-end composite wall: an insulated brick wall with inside
    /// and outside convective films, solved as a single series circuit,
    /// and cross-checked against the hand-summed resistance.
    #[test]
    fn composite_wall_end_to_end() {
        // Inside film, brick, insulation, outside film over A = 1 m^2.
        let area = 1.0;
        let r_in = convection_resistance(10.0, area).unwrap(); // 0.1
        let r_brick = conduction_resistance(0.1, 0.7, area).unwrap(); // ~0.142857
        let r_insul = conduction_resistance(0.05, 0.04, area).unwrap(); // 1.25
        let r_out = convection_resistance(40.0, area).unwrap(); // 0.025

        let net = ResistanceNetwork::series(vec![
            ResistanceNetwork::leaf(r_in).unwrap(),
            ResistanceNetwork::leaf(r_brick).unwrap(),
            ResistanceNetwork::leaf(r_insul).unwrap(),
            ResistanceNetwork::leaf(r_out).unwrap(),
        ])
        .unwrap();

        let hand_sum = r_in + r_brick + r_insul + r_out;
        assert!((net.total_resistance() - hand_sum).abs() < 1e-12);

        // Q = ΔT / R_total for a 25 °C → −5 °C drop.
        let q = net.heat_rate(25.0, -5.0).unwrap();
        assert!((q - 30.0 / hand_sum).abs() < 1e-9);
        assert!(q > 0.0);
    }

    /// Adding insulation (raising its thickness) must lower the steady
    /// heat loss through the same composite wall.
    #[test]
    fn thicker_insulation_lowers_loss() {
        let area = 1.0;
        let r_in = convection_resistance(10.0, area).unwrap();
        let r_out = convection_resistance(40.0, area).unwrap();

        let loss_with_insul = |thk: f64| -> f64 {
            let r_insul = conduction_resistance(thk, 0.04, area).unwrap();
            let net = ResistanceNetwork::series(vec![
                ResistanceNetwork::leaf(r_in).unwrap(),
                ResistanceNetwork::leaf(r_insul).unwrap(),
                ResistanceNetwork::leaf(r_out).unwrap(),
            ])
            .unwrap();
            net.heat_rate(25.0, -5.0).unwrap()
        };

        let thin = loss_with_insul(0.02);
        let thick = loss_with_insul(0.10);
        assert!(thick < thin);
    }

    /// The public re-exports line up with the underlying types.
    #[test]
    fn fin_efficiency_bounded_via_reexport() {
        let f = Fin::rectangular(0.05, 0.04, 0.002, 180.0, 25.0).unwrap();
        let eta = f.efficiency(TipCondition::Adiabatic);
        assert!(eta > 0.0 && eta <= 1.0 + EPS);
    }
}
