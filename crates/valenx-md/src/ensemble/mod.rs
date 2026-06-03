//! Thermostats and barostats — sampling the NVT and NPT ensembles.
//!
//! **Roadmap features 19–22.** A bare integrator samples the
//! microcanonical (NVE) ensemble — constant energy. Real simulations
//! want constant *temperature* (NVT) or constant temperature *and
//! pressure* (NPT). A **thermostat** couples the system to a heat
//! bath; a **barostat** couples it to a pressure bath.
//!
//! - [`berendsen`] — the Berendsen weak-coupling thermostat
//!   (feature 19).
//! - [`andersen`] — the Andersen stochastic-collision thermostat and
//!   the velocity-rescale (Bussi-Donadio-Parrinello) thermostat
//!   (feature 20).
//! - [`nose_hoover`] — the Nosé-Hoover (chain) thermostat
//!   (feature 21).
//! - [`barostat`] — the Berendsen and Parrinello-Rahman barostats
//!   (feature 22).
//!
//! Thermostats implement [`Thermostat`] and barostats [`Barostat`] so
//! the [`crate::sim`] driver can apply them generically after each
//! integration step.

pub mod andersen;
pub mod barostat;
pub mod berendsen;
pub mod nose_hoover;

use crate::error::Result;
use crate::system::System;

/// A temperature-control algorithm, applied once per step *after* the
/// integrator has advanced the system.
pub trait Thermostat {
    /// A short human-readable name for reports.
    fn name(&self) -> &str;

    /// The target temperature (K).
    fn target_temperature(&self) -> f64;

    /// Adjusts `system`'s velocities toward the target temperature.
    ///
    /// `dt` is the integration time step (ps) — coupling thermostats
    /// need it; instantaneous ones ignore it. `constraints` is the
    /// number of holonomic constraints, forwarded to the
    /// degree-of-freedom count.
    ///
    /// # Errors
    /// Implementation-specific.
    fn apply(&mut self, system: &mut System, dt: f64, constraints: usize) -> Result<()>;
}

/// A pressure-control algorithm, applied once per step.
pub trait Barostat {
    /// A short human-readable name for reports.
    fn name(&self) -> &str;

    /// The target pressure (bar).
    fn target_pressure(&self) -> f64;

    /// Rescales `system`'s box and coordinates toward the target
    /// pressure.
    ///
    /// `instantaneous_pressure` is the current virial pressure (bar)
    /// the caller has already computed; `dt` is the time step (ps).
    ///
    /// # Errors
    /// Implementation-specific.
    fn apply(
        &mut self,
        system: &mut System,
        instantaneous_pressure: f64,
        dt: f64,
    ) -> Result<()>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::berendsen::Berendsen;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    /// Shared helper: a hot ideal gas a thermostat must cool.
    pub(crate) fn hot_gas(n: usize, speed: f64) -> System {
        let mut top = Topology::new();
        for _ in 0..n {
            top.push_atom(Atom::new("A", 20.0, 0.0).unwrap());
        }
        let pos = (0..n)
            .map(|i| Vector3::new(i as f64 * 0.5, 0.0, 0.0))
            .collect();
        let mut sys = System::new(top, pos).unwrap();
        let vels = (0..n)
            .map(|i| {
                let s = if i % 2 == 0 { speed } else { -speed };
                Vector3::new(s, s * 0.5, -s * 0.3)
            })
            .collect();
        sys.set_velocities(vels).unwrap();
        sys
    }

    #[test]
    fn thermostat_trait_object_works() {
        let mut sys = hot_gas(20, 2.0);
        let mut thermo: Box<dyn Thermostat> = Box::new(Berendsen::new(300.0, 0.1).unwrap());
        thermo.apply(&mut sys, 0.002, 0).unwrap();
        assert_eq!(thermo.target_temperature(), 300.0);
    }
}
