//! Energy / temperature / pressure reporters — **roadmap feature 26**.
//!
//! The instantaneous thermodynamic state of a system, and a log that
//! records it as a run progresses:
//!
//! - **Energies** — kinetic (`Σ ½mv²`), potential (from the force
//!   evaluation) and their sum. Total energy should be conserved in an
//!   NVE run, which is the standard correctness check.
//! - **Temperature** — from equipartition, `T = 2·KE/(N_f·k_B)`.
//! - **Pressure** — from the virial theorem,
//!
//!   ```text
//!   P = [ 2·KE − Ξ ] / (3·V)
//!   ```
//!
//!   where `Ξ = −Σ rᵢⱼ·fᵢⱼ` is the (negative) inner virial that every
//!   [`crate::bonded::ForceTerm`] accumulates into
//!   [`EnergyForce::virial`]. The result is converted to bar.
//!
//! [`ObservableLog`] collects a [`StateReport`] per recorded step and
//! offers running means — what an MD report or a quick plot needs.

use crate::bonded::EnergyForce;
use crate::system::System;
use crate::units::PRESSURE_KJMOLNM3_TO_BAR;

/// A snapshot of a system's thermodynamic state.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct StateReport {
    /// Simulation step index this report was taken at.
    pub step: usize,
    /// Simulation time (ps).
    pub time: f64,
    /// Kinetic energy (kJ/mol).
    pub kinetic_energy: f64,
    /// Potential energy (kJ/mol).
    pub potential_energy: f64,
    /// Total energy (kJ/mol).
    pub total_energy: f64,
    /// Instantaneous temperature (K).
    pub temperature: f64,
    /// Instantaneous virial pressure (bar). `None` for a non-periodic
    /// system, where pressure is undefined.
    pub pressure: Option<f64>,
}

/// Computes a [`StateReport`] for a system given the potential energy
/// and accumulated virial from a force evaluation.
///
/// `constraints` is the number of holonomic constraints (forwarded to
/// the degree-of-freedom count); pass 0 if none.
pub fn state_report(
    system: &System,
    ef: &EnergyForce,
    step: usize,
    time: f64,
    constraints: usize,
) -> StateReport {
    let ke = system.kinetic_energy();
    let pe = ef.energy;
    let temperature = system.temperature(constraints);
    let pressure = pressure_bar(system, ef);
    StateReport {
        step,
        time,
        kinetic_energy: ke,
        potential_energy: pe,
        total_energy: ke + pe,
        temperature,
        pressure,
    }
}

/// The virial pressure (bar) of a periodic system, or `None` if the
/// box is non-periodic.
///
/// `P = (2·KE − Ξ)/(3V)`, with `Ξ` the inner virial. The
/// [`EnergyForce::virial`] field stores `Σ rᵢⱼ·fᵢⱼ`; for a pair force
/// that is `+Σ r·f`, and the pressure formula wants
/// `P = (2·KE + Σr·f)/(3V)` — an attractive interaction
/// (`Σr·f < 0`) lowers the pressure, as it must.
pub fn pressure_bar(system: &System, ef: &EnergyForce) -> Option<f64> {
    let volume = system.cell.volume();
    if !volume.is_finite() || volume <= 0.0 {
        return None;
    }
    let two_ke = 2.0 * system.kinetic_energy();
    let p_internal = (two_ke + ef.virial) / (3.0 * volume);
    Some(p_internal * PRESSURE_KJMOLNM3_TO_BAR)
}

/// A time series of [`StateReport`]s collected over a run.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ObservableLog {
    /// The recorded reports, in order.
    pub reports: Vec<StateReport>,
}

impl ObservableLog {
    /// An empty log.
    pub fn new() -> Self {
        ObservableLog::default()
    }

    /// Appends a report.
    pub fn record(&mut self, report: StateReport) {
        self.reports.push(report);
    }

    /// Number of recorded reports.
    pub fn len(&self) -> usize {
        self.reports.len()
    }

    /// Whether the log is empty.
    pub fn is_empty(&self) -> bool {
        self.reports.is_empty()
    }

    /// Mean total energy over the log (`0.0` if empty).
    pub fn mean_total_energy(&self) -> f64 {
        self.mean(|r| r.total_energy)
    }

    /// Mean temperature over the log.
    pub fn mean_temperature(&self) -> f64 {
        self.mean(|r| r.temperature)
    }

    /// Mean pressure over the reports that have one (`None` if no
    /// report carries a pressure).
    pub fn mean_pressure(&self) -> Option<f64> {
        let with_p: Vec<f64> = self.reports.iter().filter_map(|r| r.pressure).collect();
        if with_p.is_empty() {
            None
        } else {
            Some(with_p.iter().sum::<f64>() / with_p.len() as f64)
        }
    }

    /// Standard deviation of the total energy — a small value relative
    /// to the mean is the signature of good energy conservation.
    pub fn total_energy_std(&self) -> f64 {
        if self.reports.len() < 2 {
            return 0.0;
        }
        let mean = self.mean_total_energy();
        let var = self
            .reports
            .iter()
            .map(|r| (r.total_energy - mean).powi(2))
            .sum::<f64>()
            / self.reports.len() as f64;
        var.sqrt()
    }

    fn mean(&self, f: impl Fn(&StateReport) -> f64) -> f64 {
        if self.reports.is_empty() {
            return 0.0;
        }
        self.reports.iter().map(f).sum::<f64>() / self.reports.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::{Atom, Topology};
    use nalgebra::Vector3;

    fn moving_system() -> System {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 10.0, 0.0).unwrap());
        top.push_atom(Atom::new("B", 10.0, 0.0).unwrap());
        let mut sys = System::new(
            top,
            vec![Vector3::zeros(), Vector3::new(0.5, 0.0, 0.0)],
        )
        .unwrap()
        .with_cell(SimBox::cubic(2.0).unwrap());
        sys.set_velocities(vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        ])
        .unwrap();
        sys
    }

    #[test]
    fn state_report_sums_energies() {
        let sys = moving_system();
        let mut ef = EnergyForce::zeros(2);
        ef.energy = -5.0;
        let r = state_report(&sys, &ef, 10, 0.02, 0);
        assert!((r.kinetic_energy - sys.kinetic_energy()).abs() < 1e-12);
        assert!((r.potential_energy - (-5.0)).abs() < 1e-12);
        assert!((r.total_energy - (r.kinetic_energy - 5.0)).abs() < 1e-12);
        assert_eq!(r.step, 10);
    }

    #[test]
    fn pressure_is_none_for_open_box() {
        let mut top = Topology::new();
        top.push_atom(Atom::new("A", 1.0, 0.0).unwrap());
        let sys = System::new(top, vec![Vector3::zeros()]).unwrap();
        let ef = EnergyForce::zeros(1);
        assert!(pressure_bar(&sys, &ef).is_none());
    }

    #[test]
    fn ideal_gas_pressure_is_positive() {
        // Pure kinetic, zero virial -> positive (ideal-gas) pressure.
        let sys = moving_system();
        let ef = EnergyForce::zeros(2);
        let p = pressure_bar(&sys, &ef).unwrap();
        assert!(p > 0.0, "ideal-gas pressure = {p}");
    }

    #[test]
    fn attractive_virial_lowers_pressure() {
        let sys = moving_system();
        let mut ef = EnergyForce::zeros(2);
        // Negative virial -> attractive -> pressure drops.
        ef.virial = -100.0;
        let p_attr = pressure_bar(&sys, &ef).unwrap();
        let p_ideal = pressure_bar(&sys, &EnergyForce::zeros(2)).unwrap();
        assert!(p_attr < p_ideal);
    }

    #[test]
    fn observable_log_running_means() {
        let mut log = ObservableLog::new();
        for step in 0..10 {
            log.record(StateReport {
                step,
                time: step as f64 * 0.01,
                kinetic_energy: 5.0,
                potential_energy: -3.0,
                total_energy: 2.0,
                temperature: 300.0,
                pressure: Some(1.0),
            });
        }
        assert_eq!(log.len(), 10);
        assert!((log.mean_total_energy() - 2.0).abs() < 1e-12);
        assert!((log.mean_temperature() - 300.0).abs() < 1e-12);
        assert!((log.mean_pressure().unwrap() - 1.0).abs() < 1e-12);
        // Constant series -> zero std.
        assert!(log.total_energy_std() < 1e-12);
    }

    #[test]
    fn empty_log_is_safe() {
        let log = ObservableLog::new();
        assert!(log.is_empty());
        assert_eq!(log.mean_total_energy(), 0.0);
        assert!(log.mean_pressure().is_none());
    }
}
