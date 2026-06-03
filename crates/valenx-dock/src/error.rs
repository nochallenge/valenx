//! Docking error taxonomy.

use thiserror::Error;
use valenx_bio::format::pdbqt::PdbqtError;

/// Errors raised by [`crate::dock`] and its helpers.
#[derive(Debug, Error)]
pub enum DockError {
    /// Underlying PDBQT parse error.
    #[error("PDBQT parse: {0}")]
    Pdbqt(#[from] PdbqtError),
    /// Atom-type string in a PDBQT record wasn't a recognized AD4 code.
    #[error("atom type: {0}")]
    AtomType(String),
    /// BRANCH referenced a serial number not seen yet.
    #[error("BRANCH references unknown serial {0}")]
    UnknownBranchSerial(i32),
    /// ROOT was never opened.
    #[error("no ROOT record in ligand")]
    NoRoot,
    /// ROOT/BRANCH/END* nesting was unbalanced.
    #[error("flexibility tree imbalanced (ROOT/BRANCH not closed)")]
    FlexibilityImbalance,
    /// Receptor file was empty or had no atoms.
    #[error("receptor has no atoms")]
    EmptyReceptor,
    /// Search box edges were not strictly positive (or center had NaN).
    #[error("search box edge `{axis}` must be positive and finite, got {value}")]
    BadBox {
        /// Which axis violated the positivity / finiteness constraint
        /// ("x", "y", or "z").
        axis: &'static str,
        /// The offending edge length in Å (NaN for non-finite center
        /// components).
        value: f64,
    },
    /// Grid spacing was zero, negative, or non-finite.
    #[error("grid_spacing must be positive and finite, got {0}")]
    BadSpacing(f64),
    /// Exhaustiveness was zero or above the Vina-supported ceiling (32).
    #[error("exhaustiveness must be in 1..=32, got {0}")]
    BadExhaustiveness(u32),
    /// `num_modes` was zero — an empty output file is never useful.
    #[error("num_modes must be >= 1, got {0}")]
    BadNumModes(u32),
    /// `energy_range` was zero, negative, or non-finite.
    #[error("energy_range must be positive and finite, got {0}")]
    BadEnergyRange(f64),
    /// IO error wrapped from std::io.
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    /// Catch-all wrapper.
    #[error("dock internal error: {0}")]
    Other(String),
}

impl From<anyhow::Error> for DockError {
    fn from(e: anyhow::Error) -> Self {
        // Stringify rather than holding the anyhow::Error so that
        // `?` from anyhow-returning helpers funnels into `Other` and
        // can never accidentally hide an `io::Error` that would have
        // gone to `Io`.
        DockError::Other(e.to_string())
    }
}

/// Coarse category an LLM can use to decide who to escalate to.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// User supplied bad data (fix the input).
    Input,
    /// User-tunable knob out of range (fix the config).
    Config,
    /// Transient runtime failure (retry may help).
    Runtime,
    /// Bug in valenx-dock (file a report).
    Internal,
}

impl DockError {
    /// Stable machine-readable code. Format: `dock.<kebab_case>`.
    /// Never changes across versions; new variants get new codes.
    pub fn code(&self) -> &'static str {
        match self {
            DockError::Pdbqt(_) => "dock.pdbqt_parse",
            DockError::AtomType(_) => "dock.bad_atom_type",
            DockError::UnknownBranchSerial(_) => "dock.unknown_branch_serial",
            DockError::NoRoot => "dock.no_root",
            DockError::FlexibilityImbalance => "dock.flex_imbalance",
            DockError::EmptyReceptor => "dock.empty_receptor",
            DockError::BadBox { .. } => "dock.bad_box",
            DockError::BadSpacing(_) => "dock.bad_spacing",
            DockError::BadExhaustiveness(_) => "dock.bad_exhaustiveness",
            DockError::BadNumModes(_) => "dock.bad_num_modes",
            DockError::BadEnergyRange(_) => "dock.bad_energy_range",
            DockError::Io(_) => "dock.io",
            DockError::Other(_) => "dock.other",
        }
    }

    /// High-level classification.
    pub fn category(&self) -> ErrorCategory {
        match self {
            DockError::Pdbqt(_)
            | DockError::AtomType(_)
            | DockError::UnknownBranchSerial(_)
            | DockError::NoRoot
            | DockError::FlexibilityImbalance
            | DockError::EmptyReceptor => ErrorCategory::Input,
            DockError::BadBox { .. }
            | DockError::BadSpacing(_)
            | DockError::BadExhaustiveness(_)
            | DockError::BadNumModes(_)
            | DockError::BadEnergyRange(_) => ErrorCategory::Config,
            DockError::Io(_) => ErrorCategory::Runtime,
            DockError::Other(_) => ErrorCategory::Internal,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases = vec![
            (
                DockError::AtomType("X".into()),
                "dock.bad_atom_type",
                ErrorCategory::Input,
            ),
            (DockError::NoRoot, "dock.no_root", ErrorCategory::Input),
            (
                DockError::FlexibilityImbalance,
                "dock.flex_imbalance",
                ErrorCategory::Input,
            ),
            (
                DockError::EmptyReceptor,
                "dock.empty_receptor",
                ErrorCategory::Input,
            ),
            (
                DockError::BadBox {
                    axis: "x",
                    value: 0.0,
                },
                "dock.bad_box",
                ErrorCategory::Config,
            ),
            (
                DockError::BadSpacing(0.0),
                "dock.bad_spacing",
                ErrorCategory::Config,
            ),
            (
                DockError::BadExhaustiveness(0),
                "dock.bad_exhaustiveness",
                ErrorCategory::Config,
            ),
            (
                DockError::BadNumModes(0),
                "dock.bad_num_modes",
                ErrorCategory::Config,
            ),
            (
                DockError::BadEnergyRange(-1.0),
                "dock.bad_energy_range",
                ErrorCategory::Config,
            ),
            (
                DockError::UnknownBranchSerial(7),
                "dock.unknown_branch_serial",
                ErrorCategory::Input,
            ),
            (
                DockError::Other("oops".into()),
                "dock.other",
                ErrorCategory::Internal,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn anyhow_funnels_into_other_not_io() {
        // Important: `?` on an anyhow::Error must NOT short-circuit
        // through DockError::Io even when the underlying cause was IO.
        let e: anyhow::Error = anyhow::anyhow!("wrapped");
        let de: DockError = e.into();
        assert_eq!(de.code(), "dock.other");
    }
}
