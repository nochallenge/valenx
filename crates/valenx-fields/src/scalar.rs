//! Scalar summary records — "drag coefficient = 0.034", "max stress =
//! 120 MPa", "wall time = 42 s". Small, human-readable, queryable
//! without loading any `Field`.

use serde::{Deserialize, Serialize};

use crate::time::TimeKey;
use crate::units::Units;

/// Where the value came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ScalarSource {
    /// Computed by Valenx from raw fields (e.g. a derived
    /// coefficient).
    Computed,
    /// Extracted verbatim from the solver's output log or log-adjacent
    /// files.
    Extracted,
    /// Entered or imported from outside the solve.
    UserDefined,
}

/// One row of the results scalar catalog.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScalarRecord {
    pub name: String,
    pub value: f64,
    pub units: Units,
    pub time: TimeKey,
    pub source: ScalarSource,
    pub description: Option<String>,
}

impl ScalarRecord {
    /// Convenience constructor for the common extracted-scalar case.
    pub fn extracted(name: &str, value: f64, units: Units) -> Self {
        Self {
            name: name.to_string(),
            value,
            units,
            time: TimeKey::Steady,
            source: ScalarSource::Extracted,
            description: None,
        }
    }
}
