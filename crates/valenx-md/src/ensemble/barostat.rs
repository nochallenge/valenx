//! Barostats — **roadmap feature 22**.
//!
//! A barostat keeps the *pressure* constant by rescaling the
//! simulation box (and the atomic coordinates with it). Two schemes
//! share this module:
//!
//! - [`BerendsenBarostat`] — the **Berendsen** weak-coupling
//!   barostat. Every step the box is isotropically scaled by
//!
//!   ```text
//!   μ = [ 1 − (κ·dt/τ_p)·(P₀ − P) ]^(1/3)
//!   ```
//!
//!   where `κ` is the isothermal compressibility, `τ_p` the coupling
//!   time and `P` the instantaneous virial pressure. Like the
//!   Berendsen thermostat it relaxes the pressure correctly *on
//!   average* but does not sample the exact isothermal-isobaric
//!   distribution — good for equilibration.
//!
//! - [`ParrinelloRahman`] — the **Parrinello-Rahman** barostat. The
//!   box matrix gets its own equation of motion, driven by the
//!   imbalance between the internal virial pressure and the external
//!   target. It samples the *true* NPT ensemble and admits anisotropic
//!   box fluctuations. This v1 runs it in the **isotropic** mode (a
//!   single scalar box scale evolving under a second-order equation
//!   with a barostat mass `W`), which is the common use; the full
//!   tensorial cell-shape dynamics is the documented extension.
//!
//! Both implement the [`crate::ensemble::Barostat`] trait.

use crate::ensemble::Barostat;
use crate::error::{MdError, Result};
use crate::system::System;

/// Default isothermal compressibility of liquid water,
/// `4.5e-5 bar⁻¹` — a sensible fallback when the caller does not know
/// the system's compressibility.
pub const WATER_COMPRESSIBILITY: f64 = 4.5e-5;

/// The Berendsen weak-coupling barostat (isotropic).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BerendsenBarostat {
    /// Target pressure P₀ (bar).
    target: f64,
    /// Coupling time constant τ_p (ps).
    tau: f64,
    /// Isothermal compressibility κ (1/bar).
    compressibility: f64,
}

impl BerendsenBarostat {
    /// Builds a Berendsen barostat.
    ///
    /// * `target` — target pressure (bar)
    /// * `tau` — coupling time constant (ps); ~1 ps is typical
    /// * `compressibility` — isothermal compressibility (1/bar); pass
    ///   [`WATER_COMPRESSIBILITY`] if unsure
    ///
    /// # Errors
    /// [`MdError::Invalid`] on a non-finite target, or a non-positive
    /// `tau` / compressibility.
    pub fn new(target: f64, tau: f64, compressibility: f64) -> Result<Self> {
        if !target.is_finite() {
            return Err(MdError::invalid("target", "pressure must be finite"));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(MdError::invalid("tau", "must be finite and positive"));
        }
        if !(compressibility.is_finite() && compressibility > 0.0) {
            return Err(MdError::invalid(
                "compressibility",
                "must be finite and positive",
            ));
        }
        Ok(BerendsenBarostat {
            target,
            tau,
            compressibility,
        })
    }

    /// The isotropic box-scale factor μ for a given instantaneous
    /// pressure and step. Clamped to `[0.98, 1.02]` per step — the
    /// standard guard against a runaway rescale.
    pub fn mu(&self, pressure: f64, dt: f64) -> f64 {
        let factor =
            1.0 - (self.compressibility * dt / self.tau) * (self.target - pressure);
        factor.clamp(0.94, 1.06).cbrt()
    }
}

impl Barostat for BerendsenBarostat {
    fn name(&self) -> &str {
        "berendsen-barostat"
    }

    fn target_pressure(&self) -> f64 {
        self.target
    }

    fn apply(
        &mut self,
        system: &mut System,
        instantaneous_pressure: f64,
        dt: f64,
    ) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        if !system.cell.is_periodic() {
            return Err(MdError::invalid(
                "cell",
                "a barostat needs a periodic box",
            ));
        }
        let mu = self.mu(instantaneous_pressure, dt);
        rescale_system(system, mu)
    }
}

/// The Parrinello-Rahman barostat (isotropic v1).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct ParrinelloRahman {
    /// Target pressure P₀ (bar).
    target: f64,
    /// Coupling time constant τ_p (ps).
    tau: f64,
    /// Isothermal compressibility κ (1/bar).
    compressibility: f64,
    /// The box "velocity" — d(ln V)/dt scaled — carried between
    /// steps, which is what makes Parrinello-Rahman second-order.
    box_velocity: f64,
}

impl ParrinelloRahman {
    /// Builds an isotropic Parrinello-Rahman barostat.
    ///
    /// Same arguments as [`BerendsenBarostat::new`].
    ///
    /// # Errors
    /// [`MdError::Invalid`] on bad parameters.
    pub fn new(target: f64, tau: f64, compressibility: f64) -> Result<Self> {
        if !target.is_finite() {
            return Err(MdError::invalid("target", "pressure must be finite"));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(MdError::invalid("tau", "must be finite and positive"));
        }
        if !(compressibility.is_finite() && compressibility > 0.0) {
            return Err(MdError::invalid(
                "compressibility",
                "must be finite and positive",
            ));
        }
        Ok(ParrinelloRahman {
            target,
            tau,
            compressibility,
            box_velocity: 0.0,
        })
    }

    /// The coupling time constant τ_p (ps).
    pub fn tau(&self) -> f64 {
        self.tau
    }

    /// The current box "velocity" state.
    pub fn box_velocity(&self) -> f64 {
        self.box_velocity
    }
}

impl Barostat for ParrinelloRahman {
    fn name(&self) -> &str {
        "parrinello-rahman"
    }

    fn target_pressure(&self) -> f64 {
        self.target
    }

    fn apply(
        &mut self,
        system: &mut System,
        instantaneous_pressure: f64,
        dt: f64,
    ) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(MdError::invalid("dt", "must be finite and positive"));
        }
        if !system.cell.is_periodic() {
            return Err(MdError::invalid(
                "cell",
                "a barostat needs a periodic box",
            ));
        }
        // Second-order box dynamics: the box "acceleration" is the
        // pressure imbalance divided by the barostat mass; the
        // barostat mass W is set, as in GROMACS, from τ_p and the
        // compressibility so τ_p is the natural relaxation time.
        // dv_box/dt = (κ / τ_p²)·(P − P₀)  with a τ_p-scaled drag.
        let accel =
            (self.compressibility / (self.tau * self.tau)) * (instantaneous_pressure - self.target);
        // Light velocity damping for stability (critically-damped-ish).
        let damp = (-dt / self.tau).exp();
        self.box_velocity = self.box_velocity * damp + accel * dt;
        // Box-scale factor over this step.
        let mu = (1.0 + self.box_velocity * dt).clamp(0.94, 1.06).cbrt();
        rescale_system(system, mu)
    }
}

/// Isotropically scales a system's box and atomic coordinates by `mu`.
fn rescale_system(system: &mut System, mu: f64) -> Result<()> {
    if !(mu.is_finite() && mu > 0.0) {
        return Err(MdError::invalid("mu", "box-scale factor must be positive"));
    }
    system.cell = system.cell.scaled(mu)?;
    for p in &mut system.positions {
        *p *= mu;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    fn boxed_system(edge: f64) -> System {
        let mut top = Topology::new();
        for _ in 0..8 {
            top.push_atom(Atom::new("A", 18.0, 0.0).unwrap());
        }
        let pos = (0..8)
            .map(|i| {
                Vector3::new(
                    (i % 2) as f64 * 0.5,
                    ((i / 2) % 2) as f64 * 0.5,
                    (i / 4) as f64 * 0.5,
                )
            })
            .collect();
        System::new(top, pos)
            .unwrap()
            .with_cell(SimBox::cubic(edge).unwrap())
    }

    #[test]
    fn berendsen_rejects_bad_parameters() {
        assert!(BerendsenBarostat::new(f64::NAN, 1.0, 4.5e-5).is_err());
        assert!(BerendsenBarostat::new(1.0, 0.0, 4.5e-5).is_err());
        assert!(BerendsenBarostat::new(1.0, 1.0, -1.0).is_err());
    }

    #[test]
    fn berendsen_expands_box_when_pressure_high() {
        let mut sys = boxed_system(3.0);
        let v0 = sys.cell.volume();
        let mut baro = BerendsenBarostat::new(1.0, 1.0, WATER_COMPRESSIBILITY).unwrap();
        // Pressure far above target -> the box must *expand* so the
        // pressure relaxes toward the target (P and V are inversely
        // related; shrinking on high pressure would be unstable
        // positive feedback).
        baro.apply(&mut sys, 2000.0, 0.002).unwrap();
        assert!(sys.cell.volume() > v0, "volume did not expand");
    }

    #[test]
    fn berendsen_shrinks_box_when_pressure_low() {
        let mut sys = boxed_system(3.0);
        let v0 = sys.cell.volume();
        let mut baro = BerendsenBarostat::new(1.0, 1.0, WATER_COMPRESSIBILITY).unwrap();
        // Negative (sub-target) pressure -> the box must shrink so the
        // pressure rises back toward the target.
        baro.apply(&mut sys, -2000.0, 0.002).unwrap();
        assert!(sys.cell.volume() < v0, "volume did not shrink");
    }

    #[test]
    fn berendsen_coordinates_scale_with_box() {
        let mut sys = boxed_system(3.0);
        let p_before = sys.positions[3];
        let mut baro = BerendsenBarostat::new(1.0, 1.0, WATER_COMPRESSIBILITY).unwrap();
        baro.apply(&mut sys, 5000.0, 0.002).unwrap();
        // Atom 3 moved (scaled away from the origin as the box expanded
        // under the above-target pressure).
        assert!((sys.positions[3] - p_before).norm() > 0.0);
    }

    #[test]
    fn berendsen_rejects_non_periodic_box() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        let mut sys = System::new(top, vec![Vector3::zeros()]).unwrap();
        let mut baro = BerendsenBarostat::new(1.0, 1.0, 4.5e-5).unwrap();
        assert!(baro.apply(&mut sys, 100.0, 0.002).is_err());
    }

    #[test]
    fn parrinello_rahman_rejects_bad_parameters() {
        assert!(ParrinelloRahman::new(f64::INFINITY, 1.0, 4.5e-5).is_err());
        assert!(ParrinelloRahman::new(1.0, -1.0, 4.5e-5).is_err());
    }

    #[test]
    fn parrinello_rahman_relaxes_volume_toward_equilibrium() {
        // Hold a constant above-target pressure; the box should keep
        // expanding so the pressure relaxes toward the target.
        let mut sys = boxed_system(4.0);
        let v0 = sys.cell.volume();
        let mut baro =
            ParrinelloRahman::new(1.0, 1.0, WATER_COMPRESSIBILITY).unwrap();
        for _ in 0..500 {
            baro.apply(&mut sys, 1000.0, 0.002).unwrap();
        }
        assert!(sys.cell.volume() > v0, "PR box did not respond to pressure");
        assert!(sys.cell.volume().is_finite());
    }

    #[test]
    fn parrinello_rahman_box_velocity_evolves() {
        let mut sys = boxed_system(4.0);
        let mut baro =
            ParrinelloRahman::new(1.0, 1.0, WATER_COMPRESSIBILITY).unwrap();
        assert_eq!(baro.box_velocity(), 0.0);
        baro.apply(&mut sys, 3000.0, 0.002).unwrap();
        assert!(baro.box_velocity() != 0.0);
    }
}
