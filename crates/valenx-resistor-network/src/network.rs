//! Serializable description of a simple resistor combination.
//!
//! ## Model
//!
//! [`Combination`] is a tagged enum describing either a list of
//! resistances wired in [`Combination::Series`] or in
//! [`Combination::Parallel`]. It carries no behaviour beyond
//! [`Combination::resistance`], which dispatches to the closed forms
//! in [`crate::combination`]. The type derives
//! [`serde::Serialize`] / [`serde::Deserialize`] so a network
//! description can be round-tripped to and from a config file or
//! sent across a process boundary.

use serde::{Deserialize, Serialize};

use crate::combination::{parallel, series};
use crate::error::ResistorError;

/// A flat resistor combination: a set of resistances wired either in
/// series or in parallel.
///
/// This is deliberately shallow — it models one level of series *or*
/// parallel, which is enough to describe the building block the rest
/// of the crate computes on. Nested topologies are expressed by the
/// caller composing equivalent values (see the crate-level example).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "resistors", rename_all = "lowercase")]
pub enum Combination {
    /// Resistances wired in series; equivalent resistance is their
    /// sum.
    Series(Vec<f64>),
    /// Resistances wired in parallel; equivalent resistance is the
    /// reciprocal of the summed conductances.
    Parallel(Vec<f64>),
}

impl Combination {
    /// Equivalent resistance of this combination.
    ///
    /// Dispatches to [`series`] or [`parallel`]; see those functions
    /// for the validation rules and error cases.
    ///
    /// # Errors
    ///
    /// Propagates any [`ResistorError`] from the underlying
    /// reduction (empty network, non-positive, or non-finite arm).
    ///
    /// # Examples
    ///
    /// ```
    /// use valenx_resistor_network::network::Combination;
    /// let c = Combination::Parallel(vec![1000.0, 1000.0]);
    /// assert!((c.resistance().unwrap() - 500.0).abs() < 1e-9);
    /// ```
    pub fn resistance(&self) -> Result<f64, ResistorError> {
        match self {
            Combination::Series(rs) => series(rs),
            Combination::Parallel(rs) => parallel(rs),
        }
    }

    /// Number of resistors in the combination.
    pub fn len(&self) -> usize {
        match self {
            Combination::Series(rs) | Combination::Parallel(rs) => rs.len(),
        }
    }

    /// Whether the combination has no resistors.
    ///
    /// An empty combination is invalid input to
    /// [`Combination::resistance`]; this lets callers check before
    /// computing.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn series_combination_resistance() {
        let c = Combination::Series(vec![100.0, 220.0, 330.0]);
        assert!((c.resistance().expect("valid") - 650.0).abs() < EPS);
        assert_eq!(c.len(), 3);
        assert!(!c.is_empty());
    }

    #[test]
    fn parallel_combination_resistance() {
        let c = Combination::Parallel(vec![1000.0, 1000.0]);
        assert!((c.resistance().expect("valid") - 500.0).abs() < EPS);
    }

    #[test]
    fn empty_combination_reports_empty_and_errors() {
        let c = Combination::Series(vec![]);
        assert!(c.is_empty());
        assert_eq!(c.resistance(), Err(ResistorError::empty_network()));
    }

    #[test]
    fn round_trips_through_json() {
        let c = Combination::Parallel(vec![2000.0, 3000.0]);
        let json = serde_json::to_string(&c).expect("serialize");
        let back: Combination = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(c, back);
        assert!((back.resistance().expect("valid") - 1200.0).abs() < EPS);
    }

    #[test]
    fn serialized_shape_is_tagged() {
        let c = Combination::Series(vec![1.0, 2.0]);
        let json = serde_json::to_string(&c).expect("serialize");
        assert!(json.contains("\"kind\":\"series\""), "json was {json}");
        assert!(json.contains("\"resistors\":[1.0,2.0]"), "json was {json}");
    }
}
