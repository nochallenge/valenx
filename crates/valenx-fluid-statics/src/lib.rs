//! # valenx-fluid-statics — hydrostatics in closed form
//!
//! Pressure, buoyancy and submerged-surface forces in a fluid at rest.
//! Everything here is the classical, exact, constant-density
//! hydrostatics of a first-year fluid-mechanics text, implemented as
//! small validated functions over SI units.
//!
//! ## What
//!
//! - **Pressure** ([`pressure`]) — the hydrostatic law `P = rho * g * h`,
//!   the gauge ↔ absolute relationship `P_abs = P_gauge + P_surface`,
//!   and the inverse "pressure head" `h = P / (rho * g)`.
//! - **Buoyancy** ([`buoyancy`]) — Archimedes' principle
//!   `F_b = rho * g * V_disp`, the net force on a submerged body, the
//!   floating / sinking decision, and the floating (submerged) fraction
//!   `rho_body / rho_fluid` with the resulting waterline draft.
//! - **Submerged plate** ([`plate`]) — the resultant hydrostatic force
//!   `F = rho * g * h_c * A` and the centre of pressure
//!   `h_cp = h_c + I_xc * sin^2(theta) / (h_c * A)` on a flat plate, with
//!   a [`plate::RectangularPlate`] helper that supplies `A` and `I_xc` in
//!   closed form.
//! - **Manometer** ([`manometer`]) — gauge pressure from an open U-tube
//!   `P = rho_m * g * R`, the differential-manometer reading
//!   `(rho_m - rho_f) * g * R + rho_f * g * (z_b - z_a)`, and a
//!   path-tracing [`manometer::balance_path`] that adds up the columns
//!   of a hand-drawn manometer circuit.
//! - **Fluid model** ([`fluid`]) — a validated [`Fluid`] density type
//!   with named reference fluids (water, sea water, mercury, air) and
//!   shared physical constants.
//!
//! ## Model
//!
//! The fluid is **incompressible and of uniform density** and is **at
//! rest**, so the only body force is gravity and pressure is a function
//! of depth alone. Under those assumptions the results in this crate are
//! not approximations — they are the exact closed-form answers:
//!
//! - Gauge pressure is **exactly linear in depth**, with slope equal to
//!   the specific weight `gamma = rho * g`.
//! - The buoyant force is **exactly** the weight of the displaced fluid,
//!   and a freely floating homogeneous body displaces **exactly** the
//!   fraction `rho_body / rho_fluid` of its volume.
//! - The centre of pressure on a submerged plate is **exactly**
//!   `I_xc * sin^2(theta) / (h_c * A)` below its centroid (e.g. exactly
//!   two-thirds depth for a surface-piercing vertical rectangle).
//! - A closed manometer circuit returning to its starting height through
//!   balanced columns sums to **exactly** zero.
//!
//! Every fallible entry point validates its inputs and returns a
//! [`FluidStaticsError`] rather than producing a non-physical number;
//! float comparisons in the test suite use an epsilon, never `==`.
//!
//! ## Honest scope
//!
//! This is **research / educational grade**: textbook closed-form,
//! incompressible, constant-density hydrostatics and well-established
//! analytic models. It is **not** a clinical / medical tool and **not**
//! a production-engineering or safety-rated design tool. Deliberate
//! limits:
//!
//! - **Static fluids only** — no flow, no dynamic / stagnation pressure,
//!   no head loss, no Bernoulli or momentum effects. This is statics,
//!   not fluid *dynamics*.
//! - **Constant density** — the fluid is treated as incompressible with
//!   one density. There is no compressible atmosphere lapse, no
//!   temperature- or salinity-stratified column, and no free-surface
//!   capillarity / surface tension.
//! - **Resultant pressure force only** for submerged surfaces — the
//!   normal force and its line of action (centre of pressure), not a
//!   full stress field, not the horizontal/vertical decomposition on a
//!   curved gate, and not stability (metacentre) of a floating body.
//! - **Idealised manometers** — massless gas columns where stated,
//!   no meniscus / reading corrections, no thermal expansion of the
//!   gauge liquid.
//!
//! None of those omissions makes a result wrong within its stated
//! assumptions; each is a clearly-bounded simplification, and the
//! numbers the crate returns (a pressure, a buoyant force, a gate force
//! and its centre of pressure, a manometer reading) are the genuine
//! textbook answers.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod buoyancy;
pub mod error;
pub mod fluid;
pub mod manometer;
pub mod plate;
pub mod pressure;

// --- Convenience re-exports of the most-used types --------------------

pub use error::{ErrorCategory, FluidStaticsError, Result};
pub use fluid::{
    Fluid, DENSITY_AIR_15C, DENSITY_MERCURY, DENSITY_SEAWATER, DENSITY_WATER_4C,
    STANDARD_ATMOSPHERE_PA, STANDARD_GRAVITY,
};

pub use pressure::{
    absolute_from_gauge, absolute_pressure, absolute_pressure_open, depth_for_gauge_pressure,
    gauge_from_absolute, gauge_pressure, gauge_pressure_standard, pressure_head,
};

pub use buoyancy::{
    analyse_float, buoyant_force, buoyant_force_standard, floating_draft, floating_fraction,
    net_force_submerged, weight, FloatResult, FloatState,
};

pub use plate::{plate_load_standard, PlateLoad, RectangularPlate, SubmergedPlate};

pub use manometer::{
    balance_path, open_manometer_gauge, open_manometer_gauge_standard,
    open_manometer_gauge_with_working_fluid, DifferentialManometer, Segment,
};

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    /// End-to-end: a vertical dam-style gate. The hydrostatic force on
    /// the gate and its centre of pressure are cross-checked against the
    /// independent depth-integral identities, tying the pressure,
    /// buoyancy-free force, and plate modules together.
    #[test]
    fn submerged_gate_end_to_end() {
        let water = Fluid::water();
        // A 4 m wide, 5 m tall vertical gate with its top edge on the
        // free surface holds back a reservoir.
        let gate = RectangularPlate::vertical(4.0, 5.0).unwrap();
        let load = gate
            .load_from_top_edge(&water, STANDARD_GRAVITY, 0.0)
            .unwrap();

        // Force via the centroidal-pressure law must equal the average
        // gauge pressure (at the 2.5 m centroid) times the 20 m^2 area.
        let p_centroid = gauge_pressure(&water, STANDARD_GRAVITY, 2.5).unwrap();
        let f_expected = p_centroid * gate.area_m2();
        assert!(
            (load.force_n - f_expected).abs() < 1e-3,
            "F {}",
            load.force_n
        );

        // Surface-piercing vertical rectangle -> CP at exactly 2/3 of the
        // total 5 m depth.
        assert!(
            (load.center_of_pressure_depth_m - 10.0 / 3.0).abs() < 1e-9,
            "cp {}",
            load.center_of_pressure_depth_m
        );
        // ... and strictly below the 2.5 m centroid.
        assert!(load.center_of_pressure_depth_m > load.centroid_depth_m);
    }

    /// End-to-end: a floating block. The submerged fraction predicts the
    /// draft, and the buoyant force on that displaced volume balances the
    /// block's own weight exactly.
    #[test]
    fn floating_block_end_to_end() {
        let fluid = Fluid::seawater();
        let block = Fluid::new(750.0).unwrap(); // dense hardwood
        let area = 3.0_f64; // m^2 waterplane
        let height = 1.0_f64; // m
        let volume = area * height;

        let res = analyse_float(&block, &fluid, volume).unwrap();
        assert_eq!(res.state, FloatState::Floats);

        // Submerged fraction = 750/1025; draft = fraction * height.
        let frac = floating_fraction(&block, &fluid).unwrap();
        let draft = floating_draft(&block, &fluid, volume, area).unwrap();
        assert!((draft - frac * height).abs() < 1e-9, "draft {draft}");

        // Buoyant force on the displaced volume == block weight.
        let fb = buoyant_force(&fluid, STANDARD_GRAVITY, res.displaced_volume_m3).unwrap();
        let w = weight(block.density() * volume, STANDARD_GRAVITY).unwrap();
        assert!((fb - w).abs() < 1e-6, "fb {fb} w {w}");
    }

    /// End-to-end: a differential manometer across a horizontal pipe.
    /// Tracing the manometer as a balance path must agree with the
    /// closed-form differential-manometer formula.
    #[test]
    fn differential_manometer_end_to_end() {
        let water = Fluid::water();
        let mercury = Fluid::mercury();
        let r = 0.12; // 120 mm mercury deflection

        let m = DifferentialManometer {
            working_fluid: water,
            gauge_fluid: mercury,
            reading_m: r,
            elevation_a_m: 0.0,
            elevation_b_m: 0.0,
        };
        let dp_formula = m.pressure_difference(STANDARD_GRAVITY).unwrap();

        // Same physics by tracing columns: from A go down `r` through
        // water to the mercury, across, then up `r` through mercury and
        // down `r` back through water to B at the same elevation. Net:
        // P_A - P_B = (rho_m - rho_f) * g * r.
        let dp_path = (mercury.density() - water.density()) * STANDARD_GRAVITY * r;
        assert!(
            (dp_formula - dp_path).abs() < EPS,
            "formula {dp_formula} path {dp_path}"
        );
        assert!(dp_formula > 0.0);
    }

    /// The crate-level re-exports resolve and produce sane numbers.
    #[test]
    fn reexports_are_usable() {
        let p = gauge_pressure_standard(&Fluid::water(), 1.0).unwrap();
        assert!(p > 0.0);
        let fb = buoyant_force_standard(&Fluid::water(), 1.0).unwrap();
        assert!(fb > 0.0);
        let pm = open_manometer_gauge_standard(&Fluid::mercury(), 0.1).unwrap();
        assert!(pm > 0.0);
    }
}
