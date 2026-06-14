//! The [`Flywheel`] aggregate — a rotor plus a material density.
//!
//! [`Flywheel`] couples a validated [`Rotor`] with the material density
//! `rho` needed for the rim-stress model, and exposes the energy,
//! governor-sizing, and stress relations as methods so a caller works
//! with one object instead of threading the moment of inertia and radius
//! through free functions by hand.

use serde::{Deserialize, Serialize};

use crate::energy::{energy_fluctuation, kinetic_energy, usable_energy};
use crate::error::FlywheelError;
use crate::rotor::Rotor;
use crate::stress::{rim_speed, rim_stress};

/// A complete flywheel model: a rotor geometry plus the density of the
/// material it is made from.
///
/// Construct with [`Flywheel::new`], which validates the density. The
/// rotor is already validated by its own constructor.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Flywheel {
    /// The rotor geometry and mass.
    pub rotor: Rotor,
    /// Material density `rho` (kg/m^3), used by the rim-stress model.
    pub density: f64,
}

impl Flywheel {
    /// Build a flywheel from a [`Rotor`] and a material density (kg/m^3).
    ///
    /// # Errors
    ///
    /// Returns [`FlywheelError::InvalidParameter`] if `density` is not
    /// finite and strictly positive.
    pub fn new(rotor: Rotor, density: f64) -> Result<Self, FlywheelError> {
        let density = FlywheelError::require_positive("density", density)?;
        Ok(Self { rotor, density })
    }

    /// The rotor's mass moment of inertia about the spin axis (kg.m^2).
    #[must_use]
    pub fn moment_of_inertia(&self) -> f64 {
        self.rotor.moment_of_inertia()
    }

    /// Stored kinetic energy `E = 1/2 I omega^2` at angular speed `omega`
    /// (rad/s), in joules.
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`kinetic_energy`] for a
    /// negative / non-finite `omega`.
    pub fn energy_at(&self, omega: f64) -> Result<f64, FlywheelError> {
        kinetic_energy(self.moment_of_inertia(), omega)
    }

    /// Energy extractable as the rotor slows from `omega_max` to
    /// `omega_min`, in joules (see [`usable_energy`]).
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`usable_energy`].
    pub fn usable_energy(&self, omega_min: f64, omega_max: f64) -> Result<f64, FlywheelError> {
        usable_energy(self.moment_of_inertia(), omega_min, omega_max)
    }

    /// Energy fluctuation in the governor form `dE = I omega_avg^2 Cs`,
    /// in joules (see [`energy_fluctuation`]).
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`energy_fluctuation`].
    pub fn energy_fluctuation(&self, omega_avg: f64, cs: f64) -> Result<f64, FlywheelError> {
        energy_fluctuation(self.moment_of_inertia(), omega_avg, cs)
    }

    /// Rim (tangential) speed `v = omega r_out` at the outermost radius,
    /// in m/s (see [`rim_speed`]).
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`rim_speed`].
    pub fn rim_speed(&self, omega: f64) -> Result<f64, FlywheelError> {
        rim_speed(omega, self.rotor.outer_radius())
    }

    /// Thin-ring hoop stress `sigma = rho (omega r_out)^2` at the
    /// outermost radius, in pascals (see [`rim_stress`]).
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`rim_stress`].
    pub fn rim_stress(&self, omega: f64) -> Result<f64, FlywheelError> {
        rim_stress(self.density, omega, self.rotor.outer_radius())
    }

    /// Specific energy (energy per unit rotor mass) `E / m` at angular
    /// speed `omega`, in joules per kilogram.
    ///
    /// # Errors
    ///
    /// Propagates [`FlywheelError`] from [`Flywheel::energy_at`].
    pub fn specific_energy(&self, omega: f64) -> Result<f64, FlywheelError> {
        Ok(self.energy_at(omega)? / self.rotor.mass())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::energy::{coefficient_of_fluctuation, rpm_to_rad_s};

    const EPS: f64 = 1e-6;

    fn steel_disk() -> Flywheel {
        // 10 kg solid steel disk, 0.3 m radius.
        let rotor = Rotor::solid_disk(10.0, 0.3).unwrap();
        Flywheel::new(rotor, 7800.0).unwrap()
    }

    #[test]
    fn energy_at_matches_half_i_omega_squared() {
        let fw = steel_disk();
        // I = 0.5 * 10 * 0.09 = 0.45; omega = 100 -> E = 0.5*0.45*10000
        //                                              = 2250 J
        let e = fw.energy_at(100.0).unwrap();
        assert!((fw.moment_of_inertia() - 0.45).abs() < EPS);
        assert!((e - 2250.0).abs() < EPS);
    }

    #[test]
    fn doubling_speed_quadruples_stored_energy() {
        let fw = steel_disk();
        let e1 = fw.energy_at(120.0).unwrap();
        let e2 = fw.energy_at(240.0).unwrap();
        assert!((e2 / e1 - 4.0).abs() < 1e-9);
    }

    #[test]
    fn rim_stress_matches_rho_v_squared() {
        let fw = steel_disk();
        let omega = 100.0;
        let v = fw.rim_speed(omega).unwrap(); // 100 * 0.3 = 30 m/s
        assert!((v - 30.0).abs() < EPS);
        let sigma = fw.rim_stress(omega).unwrap();
        assert!((sigma - 7800.0 * 30.0 * 30.0).abs() < 1e-3);
    }

    #[test]
    fn governor_form_agrees_with_usable_energy() {
        // Engine flywheel: mean 1500 rpm, +/- swing giving Cs = 0.02.
        let fw = steel_disk();
        let omega_avg = rpm_to_rad_s(1500.0).unwrap();
        let cs = 0.02;
        let half = 0.5 * cs * omega_avg;
        let wmin = omega_avg - half;
        let wmax = omega_avg + half;

        // The constructed band reproduces the target Cs.
        let cs_back = coefficient_of_fluctuation(wmin, wmax).unwrap();
        assert!((cs_back - cs).abs() < 1e-12);

        let governor = fw.energy_fluctuation(omega_avg, cs).unwrap();
        let direct = fw.usable_energy(wmin, wmax).unwrap();
        assert!((governor - direct).abs() < 1e-6);
    }

    #[test]
    fn specific_energy_is_energy_over_mass() {
        let fw = steel_disk();
        let omega = 200.0;
        let se = fw.specific_energy(omega).unwrap();
        let manual = fw.energy_at(omega).unwrap() / 10.0;
        assert!((se - manual).abs() < EPS);
    }

    #[test]
    fn ring_stores_twice_the_energy_of_a_disk_same_mass_radius_speed() {
        let m = 8.0;
        let r = 0.25;
        let omega = 150.0;
        let disk = Flywheel::new(Rotor::solid_disk(m, r).unwrap(), 7800.0).unwrap();
        let ring = Flywheel::new(Rotor::thin_ring(m, r).unwrap(), 7800.0).unwrap();
        let ed = disk.energy_at(omega).unwrap();
        let er = ring.energy_at(omega).unwrap();
        assert!((er - 2.0 * ed).abs() < 1e-6);
    }

    #[test]
    fn new_rejects_bad_density() {
        let rotor = Rotor::solid_disk(1.0, 1.0).unwrap();
        assert!(Flywheel::new(rotor, 0.0).is_err());
        assert!(Flywheel::new(rotor, -1.0).is_err());
        assert!(Flywheel::new(rotor, f64::NAN).is_err());
    }

    #[test]
    fn serde_round_trips() {
        let fw = steel_disk();
        let json = serde_json::to_string(&fw).unwrap();
        let back: Flywheel = serde_json::from_str(&json).unwrap();
        assert_eq!(fw, back);
    }
}
