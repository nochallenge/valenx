//! UAS vehicle definitions and integrated performance, for **both**
//! multirotor and fixed-wing configurations.
//!
//! Each vehicle is a thin, validated assembly over the composed in-house
//! aerodynamics crates: it adds the *integration* (payload margin, range,
//! endurance, a one-call performance report) that turns a point-performance
//! crate into a design tool, while delegating the underlying physics.
//!
//! - [`MultirotorUas`] wraps [`valenx_drone::Multirotor`] (actuator-disk hover
//!   momentum theory) and optionally a [`valenx_rotor::Rotor`] blade geometry
//!   for an independent BEMT cross-check of the hover power.
//! - [`FixedWingUas`] wraps [`valenx_fixedwing::Aircraft`] (parabolic-polar
//!   point performance) and adds the electric-Breguet range.
//!
//! Battery energy is modelled as installed watt-hours `E` with a *usable
//! fraction* `f` in `(0, 1]` (depth-of-discharge × reserve), so the usable
//! energy is `E·f`.

use serde::{Deserialize, Serialize};
use valenx_drone::Multirotor;
use valenx_fixedwing::Aircraft;
use valenx_rotor::{Rotor, RotorPerformance};

use crate::error::{require_positive, require_unit_fraction, UasError};

/// Standard gravity (m/s^2). Matches the composed crates.
pub const GRAVITY: f64 = 9.806_65;
/// Nominal sea-level air density (kg/m^3, ISA 15 C).
pub const SEA_LEVEL_AIR_DENSITY: f64 = 1.225;

/// A battery, as installed energy with a usable fraction.
///
/// `usable_energy_wh = energy_wh × usable_fraction`. The usable fraction
/// captures depth-of-discharge and the landing reserve a real flight plan
/// keeps in hand; a typical value is ~0.8.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Battery {
    /// Installed (nameplate) energy `E` (watt-hours), finite and positive.
    pub energy_wh: f64,
    /// Usable fraction `f` in `(0, 1]` (depth-of-discharge × reserve).
    pub usable_fraction: f64,
}

impl Battery {
    /// Build a validated battery.
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if `energy_wh` is not finite and positive;
    /// [`UasError::OutOfUnitRange`] if `usable_fraction` is not in `(0, 1]`.
    pub fn new(energy_wh: f64, usable_fraction: f64) -> Result<Self, UasError> {
        let energy_wh = require_positive("battery energy_wh", energy_wh)?;
        let usable_fraction = require_unit_fraction("battery usable_fraction", usable_fraction)?;
        Ok(Self {
            energy_wh,
            usable_fraction,
        })
    }

    /// Usable energy `E·f` (watt-hours).
    pub fn usable_energy_wh(&self) -> f64 {
        self.energy_wh * self.usable_fraction
    }

    /// Usable energy in joules (`Wh × 3600`).
    pub fn usable_energy_j(&self) -> f64 {
        self.usable_energy_wh() * 3600.0
    }
}

// ===========================================================================
// Multirotor
// ===========================================================================

/// A multirotor UAS: a momentum-theory rotor system plus a battery, payload
/// and an electrical/propulsive efficiency for forward flight.
///
/// Hover physics is delegated to [`valenx_drone::Multirotor`]. This type adds
/// the battery, payload, the cruise efficiency, and the integration
/// (endurance, range, payload margin) that makes it a design point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MultirotorUas {
    /// The composed momentum-theory multirotor (rotor count, radius, all-up
    /// mass, figure of merit, air density).
    pub rotor: Multirotor,
    /// Battery (installed energy + usable fraction).
    pub battery: Battery,
    /// Payload mass `m_pl` (kg) carried *within* the all-up mass — used for
    /// reporting and the payload-fraction; the hover thrust is sized to the
    /// full all-up `rotor.mass_kg`.
    pub payload_kg: f64,
    /// Combined motor + ESC + propeller efficiency `η` in `(0, 1]` applied to
    /// the *forward-flight* electrical-to-useful conversion. (Hover already
    /// carries the rotor figure of merit inside [`valenx_drone::Multirotor`].)
    pub drivetrain_efficiency: f64,
}

impl MultirotorUas {
    /// Build a validated multirotor UAS.
    ///
    /// `rotor_count`, `rotor_radius_m`, `all_up_mass_kg`, `air_density` must be
    /// finite/positive; `figure_of_merit` and `drivetrain_efficiency` in
    /// `(0, 1]`; `payload_kg` finite and `>= 0` and below the all-up mass.
    ///
    /// # Errors
    ///
    /// Returns the matching [`UasError`] / wrapped [`valenx_drone::DroneError`]
    /// for any violated constraint.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        rotor_count: u32,
        rotor_radius_m: f64,
        all_up_mass_kg: f64,
        figure_of_merit: f64,
        air_density: f64,
        battery: Battery,
        payload_kg: f64,
        drivetrain_efficiency: f64,
    ) -> Result<Self, UasError> {
        // valenx-drone validates rotor_count >= 1, positive radius/mass/rho,
        // and FM in (0, 1] — reuse it rather than re-checking.
        let rotor = Multirotor::new(
            rotor_count,
            rotor_radius_m,
            all_up_mass_kg,
            figure_of_merit,
            air_density,
        )?;
        let drivetrain_efficiency =
            require_unit_fraction("drivetrain_efficiency", drivetrain_efficiency)?;
        if !(payload_kg.is_finite() && payload_kg >= 0.0) {
            return Err(UasError::NonPositive {
                quantity: "payload_kg",
                value: payload_kg,
            });
        }
        if payload_kg >= all_up_mass_kg {
            return Err(UasError::NonPositive {
                // payload cannot equal/exceed the all-up mass (no airframe left)
                quantity: "all_up_mass_kg - payload_kg",
                value: all_up_mass_kg - payload_kg,
            });
        }
        Ok(Self {
            rotor,
            battery,
            payload_kg,
            drivetrain_efficiency,
        })
    }

    /// Ideal hover power `P = T^1.5 / sqrt(2·ρ·A)` with `T = m·g` (watts),
    /// from momentum (actuator-disk) theory. Delegates to
    /// [`valenx_drone::Multirotor::ideal_hover_power`].
    pub fn ideal_hover_power_w(&self) -> f64 {
        self.rotor.ideal_hover_power()
    }

    /// Actual hover shaft power `P_ideal / FM` (watts) — the power the battery
    /// must actually supply in hover.
    pub fn hover_power_w(&self) -> f64 {
        self.rotor.actual_hover_power()
    }

    /// Hover endurance `t = E_usable / P_hover` (**seconds**).
    ///
    /// Guards the power denominator: a non-finite or non-positive hover power
    /// (impossible for a validated vehicle, but checked) yields
    /// [`UasError::NonPositive`] rather than an infinite endurance.
    pub fn hover_endurance_s(&self) -> Result<f64, UasError> {
        let p = require_positive("hover_power_w", self.hover_power_w())?;
        Ok(self.battery.usable_energy_j() / p)
    }

    /// Hover endurance in minutes.
    pub fn hover_endurance_min(&self) -> Result<f64, UasError> {
        Ok(self.hover_endurance_s()? / 60.0)
    }

    /// Forward-flight **range** (metres) at a steady cruise speed.
    ///
    /// A multirotor has no efficient lifting wing, so range is built from a
    /// steady power balance, not a Breguet wing relation. The cruise *power*
    /// is taken as the hover power scaled by a `cruise_power_factor`
    /// (induced power falls in fast forward flight but parasite power on the
    /// bluff airframe rises; the lumped factor captures the net, and `1.0`
    /// recovers "cruise at hover power"). With usable energy `E_usable` and
    /// cruise speed `V`:
    ///
    /// `endurance = (E_usable · η) / (P_hover · cruise_power_factor)`,
    /// `range = endurance · V`.
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if `cruise_speed_m_s` or `cruise_power_factor`
    /// is not finite and positive, or the resulting power is non-positive.
    pub fn cruise_range_m(
        &self,
        cruise_speed_m_s: f64,
        cruise_power_factor: f64,
    ) -> Result<f64, UasError> {
        let v = require_positive("cruise_speed_m_s", cruise_speed_m_s)?;
        let factor = require_positive("cruise_power_factor", cruise_power_factor)?;
        let cruise_power = require_positive("cruise_power_w", self.hover_power_w() * factor)?;
        let energy_to_useful = self.battery.usable_energy_j() * self.drivetrain_efficiency;
        let endurance_s = energy_to_useful / cruise_power;
        Ok(endurance_s * v)
    }

    /// Maximum **additional** payload (kg) the rotors can lift beyond the
    /// current all-up mass, from the thrust margin at a given thrust-to-weight
    /// limit.
    ///
    /// The installed thrust capability is expressed as `max_thrust_to_weight`
    /// (e.g. 2.0 means the rotors can produce twice the current weight). The
    /// extra liftable weight is `(T/W_max − 1)·W`, and the extra mass is that
    /// over `g`. A `max_thrust_to_weight <= 1` leaves no margin and returns
    /// `0.0` (never negative).
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if `max_thrust_to_weight` is not finite and
    /// positive.
    pub fn max_extra_payload_kg(&self, max_thrust_to_weight: f64) -> Result<f64, UasError> {
        let tw = require_positive("max_thrust_to_weight", max_thrust_to_weight)?;
        let margin = (tw - 1.0).max(0.0); // no negative payload
        let extra_weight_n = margin * self.rotor.weight();
        Ok(extra_weight_n / GRAVITY)
    }

    /// Independent **BEMT cross-check** of the hover shaft power, from a
    /// supplied blade geometry spun at `rpm`.
    ///
    /// Where the momentum-theory [`hover_power_w`](Self::hover_power_w) lumps
    /// all rotor losses into the figure of merit, [`valenx_rotor::Rotor`]
    /// resolves the blade element-by-element. This returns the **total**
    /// rotor-system shaft power for `rotor_count` identical rotors each
    /// producing its share of the hover thrust. It does **not** re-trim the
    /// rpm to exactly match `T = m·g`; the caller supplies an rpm and gets the
    /// power that geometry draws there, to compare against the momentum-theory
    /// estimate. Returns the per-rotor [`RotorPerformance`] too for inspection.
    ///
    /// # Errors
    ///
    /// Wraps any [`valenx_rotor::RotorError`] from the blade solve (e.g.
    /// non-convergence, bad rpm/density).
    pub fn bemt_hover_power_w(
        &self,
        blade: &Rotor,
        rpm: f64,
    ) -> Result<(f64, RotorPerformance), UasError> {
        // Hover => zero axial freestream.
        let perf = blade.solve(rpm, 0.0, self.rotor.air_density)?;
        let total_power = perf.power_w * f64::from(self.rotor.rotor_count);
        Ok((total_power, perf))
    }

    /// Compute the full multirotor performance report in one call, at a given
    /// cruise speed / power factor and thrust-to-weight limit.
    ///
    /// # Errors
    ///
    /// Propagates any error from the individual analyses above.
    pub fn performance(
        &self,
        cruise_speed_m_s: f64,
        cruise_power_factor: f64,
        max_thrust_to_weight: f64,
    ) -> Result<MultirotorPerformance, UasError> {
        Ok(MultirotorPerformance {
            all_up_mass_kg: self.rotor.mass_kg,
            disk_area_m2: self.rotor.disk_area(),
            disk_loading_pa: self.rotor.disk_loading(),
            ideal_hover_power_w: self.ideal_hover_power_w(),
            hover_power_w: self.hover_power_w(),
            hover_endurance_s: self.hover_endurance_s()?,
            cruise_speed_m_s,
            cruise_range_m: self.cruise_range_m(cruise_speed_m_s, cruise_power_factor)?,
            max_extra_payload_kg: self.max_extra_payload_kg(max_thrust_to_weight)?,
            payload_fraction: self.payload_kg / self.rotor.mass_kg,
        })
    }
}

/// A multirotor UAS's computed integrated performance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MultirotorPerformance {
    /// All-up mass `m` (kg).
    pub all_up_mass_kg: f64,
    /// Total rotor disk area `A` (m^2).
    pub disk_area_m2: f64,
    /// Disk loading `W/A` (N/m^2).
    pub disk_loading_pa: f64,
    /// Ideal hover power `T^1.5/sqrt(2·ρ·A)` (W).
    pub ideal_hover_power_w: f64,
    /// Actual hover shaft power `P_ideal/FM` (W).
    pub hover_power_w: f64,
    /// Hover endurance (s).
    pub hover_endurance_s: f64,
    /// Cruise speed the range was evaluated at (m/s).
    pub cruise_speed_m_s: f64,
    /// Forward-flight range at that cruise speed (m).
    pub cruise_range_m: f64,
    /// Maximum additional payload (kg) from the thrust margin.
    pub max_extra_payload_kg: f64,
    /// Current payload as a fraction of all-up mass.
    pub payload_fraction: f64,
}

// ===========================================================================
// Fixed-wing
// ===========================================================================

/// A fixed-wing UAS: a parabolic-polar airframe plus a battery, payload and an
/// electric drivetrain efficiency, with the electric-Breguet range.
///
/// Point performance (stall speed, drag polar, `(L/D)max`) is delegated to
/// [`valenx_fixedwing::Aircraft`]; this type adds the electric range and the
/// payload margin.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FixedWingUas {
    /// The composed point-performance airframe (wing area, mass, CLmax, AR,
    /// CD0, Oswald e, air density).
    pub aircraft: Aircraft,
    /// Battery (installed energy + usable fraction).
    pub battery: Battery,
    /// Payload mass `m_pl` (kg), `>= 0` and below the all-up mass.
    pub payload_kg: f64,
    /// Overall electric drivetrain efficiency `η` in `(0, 1]` (battery → motor
    /// → propeller useful propulsive power).
    pub drivetrain_efficiency: f64,
}

impl FixedWingUas {
    /// Build a validated fixed-wing UAS.
    ///
    /// The airframe fields are validated by [`valenx_fixedwing::Aircraft::new`]
    /// (positive wing area / mass / CLmax / AR / CD0 / density; Oswald `e` in
    /// `(0, 1]`). `drivetrain_efficiency` must be in `(0, 1]`; `payload_kg`
    /// finite, `>= 0`, and below the all-up mass.
    ///
    /// # Errors
    ///
    /// Returns the matching [`UasError`] / wrapped
    /// [`valenx_fixedwing::FixedWingError`].
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        wing_area_m2: f64,
        all_up_mass_kg: f64,
        cl_max: f64,
        aspect_ratio: f64,
        cd0: f64,
        oswald_efficiency: f64,
        air_density: f64,
        battery: Battery,
        payload_kg: f64,
        drivetrain_efficiency: f64,
    ) -> Result<Self, UasError> {
        let aircraft = Aircraft::new(
            wing_area_m2,
            all_up_mass_kg,
            cl_max,
            aspect_ratio,
            cd0,
            oswald_efficiency,
            air_density,
        )?;
        let drivetrain_efficiency =
            require_unit_fraction("drivetrain_efficiency", drivetrain_efficiency)?;
        if !(payload_kg.is_finite() && payload_kg >= 0.0) {
            return Err(UasError::NonPositive {
                quantity: "payload_kg",
                value: payload_kg,
            });
        }
        if payload_kg >= all_up_mass_kg {
            return Err(UasError::NonPositive {
                quantity: "all_up_mass_kg - payload_kg",
                value: all_up_mass_kg - payload_kg,
            });
        }
        Ok(Self {
            aircraft,
            battery,
            payload_kg,
            drivetrain_efficiency,
        })
    }

    /// Maximum lift-to-drag ratio `(L/D)max`, delegated to
    /// [`valenx_fixedwing::Aircraft::max_lift_to_drag`].
    pub fn max_lift_to_drag(&self) -> f64 {
        self.aircraft.max_lift_to_drag()
    }

    /// Stall speed (m/s), delegated to
    /// [`valenx_fixedwing::Aircraft::stall_speed`].
    pub fn stall_speed_m_s(&self) -> f64 {
        self.aircraft.stall_speed()
    }

    /// **Electric-Breguet** range (metres):
    ///
    /// `R = (E_usable / (g · m)) · (L/D) · η`
    ///
    /// with `E_usable` the usable battery energy (J), `m` the all-up mass,
    /// `L/D` the lift-to-drag ratio flown, and `η` the overall electric
    /// drivetrain efficiency. Unlike the fuel-burn Breguet relation, an
    /// electric aircraft's mass is constant, so the range is *linear* in the
    /// stored energy and reduces to this closed form. If `lift_to_drag` is not
    /// supplied (`None`), the best `(L/D)max` is used.
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if the supplied / computed `L/D` or the
    /// all-up mass is not finite and positive.
    pub fn breguet_range_m(&self, lift_to_drag: Option<f64>) -> Result<f64, UasError> {
        let ld = require_positive(
            "lift_to_drag",
            lift_to_drag.unwrap_or_else(|| self.max_lift_to_drag()),
        )?;
        let mass = require_positive("all_up_mass_kg", self.aircraft.mass_kg)?;
        let energy_j = self.battery.usable_energy_j();
        // R = (E / (g m)) * (L/D) * eta.  g m > 0 guaranteed by the guards.
        Ok((energy_j / (GRAVITY * mass)) * ld * self.drivetrain_efficiency)
    }

    /// Endurance (seconds) at a steady cruise speed flown at a given `L/D`.
    ///
    /// From the same energy balance: the propulsive power needed to overcome
    /// drag `D = W/(L/D)` at speed `V` is `P = D·V/η`, so
    /// `endurance = E_usable / P = E_usable·η / (W·V/(L/D))`.
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if `cruise_speed_m_s` or the `L/D` is not
    /// finite and positive.
    pub fn endurance_s(
        &self,
        cruise_speed_m_s: f64,
        lift_to_drag: Option<f64>,
    ) -> Result<f64, UasError> {
        let v = require_positive("cruise_speed_m_s", cruise_speed_m_s)?;
        let ld = require_positive(
            "lift_to_drag",
            lift_to_drag.unwrap_or_else(|| self.max_lift_to_drag()),
        )?;
        let weight = self.aircraft.weight();
        let drag = weight / ld; // N
        let power = require_positive("cruise_power_w", drag * v / self.drivetrain_efficiency)?;
        Ok(self.battery.usable_energy_j() / power)
    }

    /// Maximum **additional** payload (kg) before the wing can no longer
    /// generate the lift to hold level flight at a chosen cruise speed.
    ///
    /// At speed `V` the wing can lift at most `L_max = 0.5·ρ·V²·S·CLmax`; the
    /// margin over the current weight, divided by `g`, is the extra mass:
    /// `(L_max − W)/g`, floored at `0.0`. `V` must be at or above the current
    /// stall speed for this to be meaningful (a slower `V` simply yields a
    /// small or zero margin, never negative).
    ///
    /// # Errors
    ///
    /// [`UasError::NonPositive`] if `cruise_speed_m_s` is not finite and
    /// positive.
    pub fn max_extra_payload_kg(&self, cruise_speed_m_s: f64) -> Result<f64, UasError> {
        let v = require_positive("cruise_speed_m_s", cruise_speed_m_s)?;
        let l_max = self.aircraft.lift(v, self.aircraft.cl_max); // N
        let margin_n = (l_max - self.aircraft.weight()).max(0.0);
        Ok(margin_n / GRAVITY)
    }

    /// Compute the full fixed-wing performance report in one call.
    ///
    /// # Errors
    ///
    /// Propagates any error from the individual analyses above.
    pub fn performance(&self, cruise_speed_m_s: f64) -> Result<FixedWingPerformance, UasError> {
        Ok(FixedWingPerformance {
            all_up_mass_kg: self.aircraft.mass_kg,
            wing_loading_pa: self.aircraft.wing_loading(),
            stall_speed_m_s: self.stall_speed_m_s(),
            max_lift_to_drag: self.max_lift_to_drag(),
            breguet_range_m: self.breguet_range_m(None)?,
            cruise_speed_m_s,
            endurance_s: self.endurance_s(cruise_speed_m_s, None)?,
            max_extra_payload_kg: self.max_extra_payload_kg(cruise_speed_m_s)?,
            payload_fraction: self.payload_kg / self.aircraft.mass_kg,
        })
    }
}

/// A fixed-wing UAS's computed integrated performance.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FixedWingPerformance {
    /// All-up mass `m` (kg).
    pub all_up_mass_kg: f64,
    /// Wing loading `W/S` (N/m^2).
    pub wing_loading_pa: f64,
    /// Stall speed `Vs` (m/s).
    pub stall_speed_m_s: f64,
    /// Maximum lift-to-drag ratio.
    pub max_lift_to_drag: f64,
    /// Electric-Breguet range at best `(L/D)` (m).
    pub breguet_range_m: f64,
    /// Cruise speed the endurance was evaluated at (m/s).
    pub cruise_speed_m_s: f64,
    /// Endurance at that cruise speed (s).
    pub endurance_s: f64,
    /// Maximum additional payload (kg) at that cruise speed.
    pub max_extra_payload_kg: f64,
    /// Current payload as a fraction of all-up mass.
    pub payload_fraction: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use valenx_rotor::Rotor;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-6 * b.abs().max(1.0)
    }

    fn batt() -> Battery {
        // 100 Wh, 80% usable.
        Battery::new(100.0, 0.8).unwrap()
    }

    /// A 1.5 kg quadcopter with 0.15 m rotors at sea level.
    fn quad() -> MultirotorUas {
        MultirotorUas::new(4, 0.15, 1.5, 0.7, SEA_LEVEL_AIR_DENSITY, batt(), 0.3, 0.75).unwrap()
    }

    /// A small fixed-wing UAS.
    fn wing() -> FixedWingUas {
        FixedWingUas::new(
            0.5,                   // wing area m^2
            3.0,                   // all-up mass kg
            1.3,                   // CLmax
            10.0,                  // aspect ratio
            0.03,                  // CD0
            0.85,                  // Oswald e
            SEA_LEVEL_AIR_DENSITY, // rho
            batt(),
            0.5,  // payload kg
            0.65, // drivetrain efficiency
        )
        .unwrap()
    }

    // ---- BENCHMARK PIN: hover power = T^1.5 / sqrt(2 rho A) ----------------
    #[test]
    fn hover_power_matches_momentum_theory_closed_form() {
        let q = quad();
        let t = q.rotor.mass_kg * GRAVITY; // T = m g
        let a = q.rotor.disk_area();
        let expected_ideal = t.powf(1.5) / (2.0 * SEA_LEVEL_AIR_DENSITY * a).sqrt();
        assert!(
            close(q.ideal_hover_power_w(), expected_ideal),
            "ideal hover power {} != closed form {}",
            q.ideal_hover_power_w(),
            expected_ideal
        );
        // Actual = ideal / FM.
        assert!(close(q.hover_power_w(), expected_ideal / 0.7));
    }

    // ---- BENCHMARK PIN: endurance = E / P ---------------------------------
    #[test]
    fn hover_endurance_is_usable_energy_over_power() {
        let q = quad();
        let expected = q.battery.usable_energy_j() / q.hover_power_w();
        assert!(close(q.hover_endurance_s().unwrap(), expected));
        // Twice the usable energy => twice the endurance.
        let q2 = MultirotorUas::new(
            4,
            0.15,
            1.5,
            0.7,
            SEA_LEVEL_AIR_DENSITY,
            Battery::new(200.0, 0.8).unwrap(),
            0.3,
            0.75,
        )
        .unwrap();
        assert!(close(
            q2.hover_endurance_s().unwrap(),
            2.0 * q.hover_endurance_s().unwrap()
        ));
    }

    #[test]
    fn endurance_matches_drone_crate_minutes() {
        // Cross-check our seconds-based endurance against valenx-drone's own
        // minutes helper on the usable energy — they must agree.
        let q = quad();
        let drone_min = q
            .rotor
            .hover_endurance_minutes(q.battery.usable_energy_wh());
        assert!(close(q.hover_endurance_min().unwrap(), drone_min));
    }

    #[test]
    fn bemt_cross_check_is_same_order_as_momentum_theory() {
        // A blade geometry roughly sized for the quad's rotor radius; the BEMT
        // total power at a hover-ish rpm should land in the same ballpark as
        // the momentum-theory estimate (same physical rotor, two methods),
        // not differ by orders of magnitude.
        let r_tip = 0.15;
        let r_hub = 0.02;
        let radii = [0.03, 0.06, 0.09, 0.12, 0.15];
        let chords = [0.025, 0.022, 0.018, 0.014, 0.008];
        let twists = [
            22.0_f64.to_radians(),
            16.0_f64.to_radians(),
            12.0_f64.to_radians(),
            9.0_f64.to_radians(),
            7.0_f64.to_radians(),
        ];
        let blade = Rotor::from_slices(2, r_tip, r_hub, &radii, &chords, &twists).unwrap();
        let q = quad();
        let (total_power, perf) = q.bemt_hover_power_w(&blade, 6500.0).unwrap();
        assert!(total_power.is_finite() && total_power > 0.0);
        assert!(perf.thrust_n.is_finite());
        // Same order of magnitude as the momentum-theory hover power.
        let mt = q.hover_power_w();
        let ratio = total_power / mt;
        assert!(
            (0.05..20.0).contains(&ratio),
            "BEMT total power {total_power} W vs momentum theory {mt} W: ratio {ratio} not within a decade-ish band"
        );
    }

    #[test]
    fn multirotor_cruise_range_is_positive_and_scales_with_speed() {
        let q = quad();
        let r1 = q.cruise_range_m(10.0, 1.0).unwrap();
        let r2 = q.cruise_range_m(20.0, 1.0).unwrap();
        assert!(r1 > 0.0 && r2 > 0.0);
        // At the same power factor, range is linear in speed.
        assert!(close(r2, 2.0 * r1));
    }

    #[test]
    fn multirotor_payload_margin_from_thrust_to_weight() {
        let q = quad();
        // T/W = 2 => can lift one extra current-weight worth of mass.
        let extra = q.max_extra_payload_kg(2.0).unwrap();
        assert!(close(extra, q.rotor.mass_kg)); // (2-1)*W/g = m
                                                // No margin at or below T/W = 1.
        assert_eq!(q.max_extra_payload_kg(1.0).unwrap(), 0.0);
        assert_eq!(q.max_extra_payload_kg(0.5).unwrap(), 0.0);
    }

    // ---- BENCHMARK PIN: fixed-wing electric-Breguet R = (E/(g m))*(L/D)*eta
    #[test]
    fn fixed_wing_breguet_matches_closed_form() {
        let w = wing();
        let ld = w.max_lift_to_drag();
        let expected = (w.battery.usable_energy_j() / (GRAVITY * w.aircraft.mass_kg)) * ld * 0.65;
        assert!(close(w.breguet_range_m(None).unwrap(), expected));
        // Explicit L/D path agrees with a hand value.
        let expected_10 =
            (w.battery.usable_energy_j() / (GRAVITY * w.aircraft.mass_kg)) * 10.0 * 0.65;
        assert!(close(w.breguet_range_m(Some(10.0)).unwrap(), expected_10));
    }

    #[test]
    fn fixed_wing_range_beats_multirotor_for_same_battery() {
        // The whole point of a wing: far more range per Wh than a hovering
        // multirotor. Sanity, not a tight number.
        let w = wing();
        let q = quad();
        let wing_range = w.breguet_range_m(None).unwrap();
        let multi_range = q.cruise_range_m(15.0, 1.0).unwrap();
        assert!(
            wing_range > multi_range,
            "fixed-wing range {wing_range} should exceed multirotor {multi_range}"
        );
    }

    #[test]
    fn fixed_wing_endurance_positive_and_falls_with_speed() {
        let w = wing();
        let vs = w.stall_speed_m_s();
        let e_slow = w.endurance_s(vs * 1.3, None).unwrap();
        let e_fast = w.endurance_s(vs * 2.5, None).unwrap();
        assert!(e_slow > 0.0 && e_fast > 0.0);
        // At fixed L/D, power ∝ V, so endurance falls as speed rises.
        assert!(e_fast < e_slow);
    }

    #[test]
    fn fixed_wing_payload_margin_grows_with_speed() {
        let w = wing();
        let vs = w.stall_speed_m_s();
        // Just above stall: little or no margin.
        let near_stall = w.max_extra_payload_kg(vs * 1.01).unwrap();
        // Much faster: a real positive margin.
        let fast = w.max_extra_payload_kg(vs * 2.0).unwrap();
        assert!(fast > near_stall);
        assert!(fast > 0.0);
        assert!(near_stall >= 0.0); // never negative
    }

    #[test]
    fn reports_round_trip_through_serde() {
        let qp = quad().performance(15.0, 1.0, 2.0).unwrap();
        let json = serde_json::to_string(&qp).unwrap();
        let back: MultirotorPerformance = serde_json::from_str(&json).unwrap();
        assert!(close(back.hover_power_w, qp.hover_power_w));
        assert!(close(back.cruise_range_m, qp.cruise_range_m));

        let wp = wing().performance(18.0).unwrap();
        let json = serde_json::to_string(&wp).unwrap();
        let back: FixedWingPerformance = serde_json::from_str(&json).unwrap();
        assert!(close(back.breguet_range_m, wp.breguet_range_m));
        assert!(close(back.max_lift_to_drag, wp.max_lift_to_drag));
    }

    // ---- Fail-loud: degenerate inputs -> Err, never a panic ---------------
    #[test]
    fn rejects_out_of_domain_construction() {
        // zero / negative rotor radius, zero mass, FM out of range — handled
        // by valenx-drone and wrapped.
        assert!(MultirotorUas::new(4, 0.0, 1.5, 0.7, 1.225, batt(), 0.3, 0.75).is_err());
        assert!(MultirotorUas::new(0, 0.15, 1.5, 0.7, 1.225, batt(), 0.3, 0.75).is_err());
        assert!(MultirotorUas::new(4, 0.15, 0.0, 0.7, 1.225, batt(), 0.3, 0.75).is_err());
        assert!(MultirotorUas::new(4, 0.15, 1.5, 1.5, 1.225, batt(), 0.3, 0.75).is_err());
        // bad drivetrain efficiency / payload >= mass / non-finite payload.
        assert!(MultirotorUas::new(4, 0.15, 1.5, 0.7, 1.225, batt(), 0.3, 0.0).is_err());
        assert!(MultirotorUas::new(4, 0.15, 1.5, 0.7, 1.225, batt(), 1.5, 0.75).is_err());
        assert!(MultirotorUas::new(4, 0.15, 1.5, 0.7, 1.225, batt(), f64::NAN, 0.75).is_err());
        // bad battery.
        assert!(Battery::new(0.0, 0.8).is_err());
        assert!(Battery::new(100.0, 0.0).is_err());
        assert!(Battery::new(100.0, 1.5).is_err());
        assert!(Battery::new(f64::INFINITY, 0.8).is_err());

        // Fixed-wing bad inputs.
        assert!(
            FixedWingUas::new(0.0, 3.0, 1.3, 10.0, 0.03, 0.85, 1.225, batt(), 0.5, 0.65).is_err()
        );
        assert!(
            FixedWingUas::new(0.5, 3.0, 1.3, 10.0, 0.03, 1.5, 1.225, batt(), 0.5, 0.65).is_err()
        );
        assert!(
            FixedWingUas::new(0.5, 3.0, 1.3, 10.0, 0.03, 0.85, 1.225, batt(), 5.0, 0.65).is_err()
        );
        assert!(
            FixedWingUas::new(0.5, 3.0, 1.3, 10.0, 0.03, 0.85, 1.225, batt(), 0.5, 0.0).is_err()
        );
    }

    #[test]
    fn rejects_degenerate_analysis_inputs() {
        let q = quad();
        assert!(q.cruise_range_m(0.0, 1.0).is_err()); // zero speed
        assert!(q.cruise_range_m(-5.0, 1.0).is_err()); // negative speed
        assert!(q.cruise_range_m(10.0, 0.0).is_err()); // zero power factor
        assert!(q.cruise_range_m(f64::NAN, 1.0).is_err()); // NaN
        assert!(q.max_extra_payload_kg(0.0).is_err()); // zero T/W
        assert!(q.max_extra_payload_kg(f64::INFINITY).is_err()); // inf T/W

        let w = wing();
        assert!(w.breguet_range_m(Some(0.0)).is_err()); // zero L/D
        assert!(w.breguet_range_m(Some(-3.0)).is_err()); // negative L/D
        assert!(w.endurance_s(0.0, None).is_err()); // zero speed
        assert!(w.max_extra_payload_kg(0.0).is_err()); // zero speed
    }
}
