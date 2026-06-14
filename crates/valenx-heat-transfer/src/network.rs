//! Thermal-resistance networks (the electrical analogue).
//!
//! ## Model
//!
//! A 1D heat path can be modelled as a circuit of thermal resistors,
//! exactly like Ohm's law with `Q` (heat rate) in place of current and
//! `ΔT` (temperature difference) in place of voltage:
//!
//! ```text
//! Q = ΔT / R_total
//! ```
//!
//! Resistors combine by the usual rules:
//!
//! - **Series** (layers stacked along the heat path): resistances add,
//!   `R_series = Σ Rᵢ`.
//! - **Parallel** (alternative paths between the same two nodes):
//!   conductances add, `1 / R_parallel = Σ (1 / Rᵢ)`.
//!
//! A [`ResistanceNetwork`] is a tree of these two combinators plus leaf
//! resistors, so an arbitrary composite wall (e.g. brick + insulation
//! in series, with a stud bridging the insulation in parallel) reduces
//! to a single equivalent resistance and hence a single heat rate.
//!
//! These reduction rules are standard (Incropera §3.1.4) and hold for
//! steady 1D conduction/convection where each branch sees the same end
//! temperatures.

use serde::{Deserialize, Serialize};

use crate::error::{require_finite, require_positive, HeatTransferError, Result};

/// A thermal-resistance network node.
///
/// Build leaves with [`ResistanceNetwork::leaf`] and combine them with
/// [`ResistanceNetwork::series`] / [`ResistanceNetwork::parallel`].
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum ResistanceNetwork {
    /// A single resistor of the given value (K/W).
    Resistor(f64),
    /// Resistors in series — total resistance is the sum.
    Series(Vec<ResistanceNetwork>),
    /// Resistors in parallel — total conductance is the sum.
    Parallel(Vec<ResistanceNetwork>),
}

impl ResistanceNetwork {
    /// A single leaf resistor of `r` K/W.
    ///
    /// # Errors
    ///
    /// Returns an error if `r` is not finite and strictly positive.
    pub fn leaf(r: f64) -> Result<Self> {
        Ok(Self::Resistor(require_positive("resistance", r)?))
    }

    /// Combine the given branches in series.
    ///
    /// # Errors
    ///
    /// Returns [`HeatTransferError::EmptyNetwork`] if `branches` is
    /// empty.
    pub fn series(branches: Vec<ResistanceNetwork>) -> Result<Self> {
        if branches.is_empty() {
            return Err(HeatTransferError::EmptyNetwork("series"));
        }
        Ok(Self::Series(branches))
    }

    /// Combine the given branches in parallel.
    ///
    /// # Errors
    ///
    /// Returns [`HeatTransferError::EmptyNetwork`] if `branches` is
    /// empty.
    pub fn parallel(branches: Vec<ResistanceNetwork>) -> Result<Self> {
        if branches.is_empty() {
            return Err(HeatTransferError::EmptyNetwork("parallel"));
        }
        Ok(Self::Parallel(branches))
    }

    /// Reduce the whole network to a single equivalent resistance
    /// (K/W).
    ///
    /// Series branches add their resistances; parallel branches add
    /// their conductances (reciprocals) and the total is the reciprocal
    /// of that sum.
    pub fn total_resistance(&self) -> f64 {
        match self {
            ResistanceNetwork::Resistor(r) => *r,
            ResistanceNetwork::Series(branches) => {
                branches.iter().map(|b| b.total_resistance()).sum()
            }
            ResistanceNetwork::Parallel(branches) => {
                let sum_g: f64 = branches.iter().map(|b| 1.0 / b.total_resistance()).sum();
                1.0 / sum_g
            }
        }
    }

    /// Steady heat rate `Q = ΔT / R_total` (W) for a temperature
    /// difference imposed across the whole network.
    ///
    /// # Errors
    ///
    /// Returns an error if either end temperature is non-finite.
    pub fn heat_rate(&self, t_hot: f64, t_cold: f64) -> Result<f64> {
        let t_hot = require_finite("t_hot", t_hot)?;
        let t_cold = require_finite("t_cold", t_cold)?;
        Ok((t_hot - t_cold) / self.total_resistance())
    }
}

/// Total resistance of resistors in series: `R = Σ Rᵢ`.
///
/// # Errors
///
/// Returns an error if `resistances` is empty or any value is not
/// finite and strictly positive.
pub fn series_resistance(resistances: &[f64]) -> Result<f64> {
    if resistances.is_empty() {
        return Err(HeatTransferError::EmptyNetwork("series"));
    }
    let mut sum = 0.0;
    for &r in resistances {
        sum += require_positive("resistance", r)?;
    }
    Ok(sum)
}

/// Equivalent resistance of resistors in parallel:
/// `1 / R = Σ (1 / Rᵢ)`.
///
/// # Errors
///
/// Returns an error if `resistances` is empty or any value is not
/// finite and strictly positive.
pub fn parallel_resistance(resistances: &[f64]) -> Result<f64> {
    if resistances.is_empty() {
        return Err(HeatTransferError::EmptyNetwork("parallel"));
    }
    let mut sum_g = 0.0;
    for &r in resistances {
        sum_g += 1.0 / require_positive("resistance", r)?;
    }
    Ok(1.0 / sum_g)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    #[test]
    fn series_resistances_add() {
        // 1 + 2 + 3 = 6 K/W.
        assert!((series_resistance(&[1.0, 2.0, 3.0]).unwrap() - 6.0).abs() < EPS);
    }

    #[test]
    fn parallel_two_equal_halves() {
        // Two 4 K/W resistors in parallel -> 2 K/W.
        assert!((parallel_resistance(&[4.0, 4.0]).unwrap() - 2.0).abs() < EPS);
    }

    #[test]
    fn parallel_is_less_than_smallest_branch() {
        let rp = parallel_resistance(&[2.0, 3.0, 6.0]).unwrap();
        // 1/R = 1/2 + 1/3 + 1/6 = 1 -> R = 1 K/W < min branch (2).
        assert!((rp - 1.0).abs() < EPS);
        assert!(rp < 2.0);
    }

    #[test]
    fn network_tree_reduces_correctly() {
        // Series( 1, Parallel(4,4), 2 ) = 1 + 2 + 2 = 5 K/W.
        let net = ResistanceNetwork::series(vec![
            ResistanceNetwork::leaf(1.0).unwrap(),
            ResistanceNetwork::parallel(vec![
                ResistanceNetwork::leaf(4.0).unwrap(),
                ResistanceNetwork::leaf(4.0).unwrap(),
            ])
            .unwrap(),
            ResistanceNetwork::leaf(2.0).unwrap(),
        ])
        .unwrap();
        assert!((net.total_resistance() - 5.0).abs() < EPS);
    }

    #[test]
    fn heat_rate_uses_total_resistance() {
        let net = ResistanceNetwork::series(vec![
            ResistanceNetwork::leaf(1.0).unwrap(),
            ResistanceNetwork::leaf(3.0).unwrap(),
        ])
        .unwrap();
        // R_total = 4, ΔT = 80 -> Q = 20 W.
        let q = net.heat_rate(100.0, 20.0).unwrap();
        assert!((q - 20.0).abs() < 1e-9);
    }

    #[test]
    fn series_then_heat_rate_matches_manual() {
        // A composite wall: conduction 0.0125 + convection 0.02 in series.
        let r_total = series_resistance(&[0.0125, 0.02]).unwrap();
        assert!((r_total - 0.0325).abs() < EPS);
        let q = (50.0 - 0.0) / r_total;
        let net = ResistanceNetwork::series(vec![
            ResistanceNetwork::leaf(0.0125).unwrap(),
            ResistanceNetwork::leaf(0.02).unwrap(),
        ])
        .unwrap();
        assert!((net.heat_rate(50.0, 0.0).unwrap() - q).abs() < 1e-9);
    }

    #[test]
    fn empty_networks_are_rejected() {
        assert!(series_resistance(&[]).is_err());
        assert!(parallel_resistance(&[]).is_err());
        assert!(ResistanceNetwork::series(vec![]).is_err());
        assert!(ResistanceNetwork::parallel(vec![]).is_err());
    }

    #[test]
    fn non_positive_leaf_rejected() {
        assert!(ResistanceNetwork::leaf(0.0).is_err());
        assert!(series_resistance(&[1.0, -2.0]).is_err());
    }
}
