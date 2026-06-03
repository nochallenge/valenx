//! Berendsen weak-coupling thermostat — **roadmap feature 19**.
//!
//! The Berendsen thermostat rescales the velocities every step by a
//! factor `λ` that nudges the temperature toward the target with a
//! first-order relaxation:
//!
//! ```text
//! dT/dt = (T₀ − T) / τ
//! λ = √( 1 + (dt/τ)·(T₀/T − 1) )
//! ```
//!
//! `τ` is the coupling time constant (ps): large `τ` couples weakly
//! (slow relaxation, minimal perturbation), small `τ` couples tightly.
//!
//! ## Honest caveat
//!
//! Berendsen relaxes the temperature *correctly on average* but it
//! does **not** sample a true canonical ensemble — it suppresses the
//! physical fluctuations of the kinetic energy ("flying ice cube" /
//! wrong `⟨ΔKE²⟩`). It is excellent for **equilibration** and is
//! provided for that; for production canonical sampling prefer the
//! velocity-rescale ([`crate::ensemble::andersen`]) or Nosé-Hoover
//! ([`crate::ensemble::nose_hoover`]) thermostats. This limitation is
//! the Berendsen scheme's, not an implementation shortcut.

use crate::ensemble::Thermostat;
use crate::error::{MdError, Result};
use crate::system::System;

/// The Berendsen weak-coupling thermostat.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Berendsen {
    /// Target temperature T₀ (K).
    target: f64,
    /// Coupling time constant τ (ps).
    tau: f64,
}

impl Berendsen {
    /// Builds a Berendsen thermostat.
    ///
    /// * `target` — target temperature (K)
    /// * `tau` — coupling time constant (ps); 0.1–1 ps is typical
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite / negative temperature or
    /// a non-positive `tau`.
    pub fn new(target: f64, tau: f64) -> Result<Self> {
        if !(target.is_finite() && target >= 0.0) {
            return Err(MdError::invalid(
                "target",
                "temperature must be finite and non-negative",
            ));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(MdError::invalid("tau", "must be finite and positive"));
        }
        Ok(Berendsen { target, tau })
    }

    /// The coupling time constant τ (ps).
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// The velocity-rescale factor λ for a current temperature `t` and
    /// step `dt`. Exposed for tests; clamped so a single step cannot
    /// rescale velocities by more than a factor of 2 (the standard
    /// stability guard).
    pub fn lambda(&self, t: f64, dt: f64) -> f64 {
        if t <= 1e-12 {
            // No kinetic energy to rescale; leave it (a noise source,
            // if any, will inject some next step).
            return 1.0;
        }
        let factor = 1.0 + (dt / self.tau) * (self.target / t - 1.0);
        factor.clamp(0.25, 4.0).sqrt()
    }
}

impl Thermostat for Berendsen {
    fn name(&self) -> &str {
        "berendsen"
    }

    fn target_temperature(&self) -> f64 {
        self.target
    }

    fn apply(&mut self, system: &mut System, dt: f64, constraints: usize) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        let t = system.temperature(constraints);
        let lambda = self.lambda(t, dt);
        for v in &mut system.velocities {
            *v *= lambda;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ensemble::tests::hot_gas;

    #[test]
    fn rejects_bad_parameters() {
        assert!(Berendsen::new(-1.0, 0.1).is_err());
        assert!(Berendsen::new(300.0, 0.0).is_err());
        assert!(Berendsen::new(300.0, -0.1).is_err());
    }

    #[test]
    fn hot_system_cools_toward_target() {
        let mut sys = hot_gas(100, 3.0);
        let t_start = sys.temperature(0);
        let mut thermo = Berendsen::new(300.0, 0.1).unwrap();
        // Should be hot to begin with.
        assert!(t_start > 300.0);
        for _ in 0..2000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
        }
        let t_end = sys.temperature(0);
        assert!((t_end - 300.0).abs() < 5.0, "T ended at {t_end}");
    }

    #[test]
    fn cold_system_heats_toward_target() {
        let mut sys = hot_gas(100, 0.2); // cold
        let t_start = sys.temperature(0);
        assert!(t_start < 300.0);
        let mut thermo = Berendsen::new(300.0, 0.1).unwrap();
        for _ in 0..3000 {
            thermo.apply(&mut sys, 0.002, 0).unwrap();
        }
        assert!((sys.temperature(0) - 300.0).abs() < 5.0);
    }

    #[test]
    fn lambda_is_one_at_target_temperature() {
        let thermo = Berendsen::new(300.0, 0.5).unwrap();
        assert!((thermo.lambda(300.0, 0.002) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn zero_temperature_is_handled() {
        let mut sys = hot_gas(10, 0.0); // exactly zero velocity
        let mut thermo = Berendsen::new(300.0, 0.1).unwrap();
        // Must not panic / NaN.
        thermo.apply(&mut sys, 0.002, 0).unwrap();
        assert!(sys.temperature(0).is_finite());
    }
}
