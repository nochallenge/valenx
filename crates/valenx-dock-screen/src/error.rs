//! Error taxonomy for `valenx-dock-screen`.
//!
//! Every fallible public function in this crate returns
//! [`Result<_, DockScreenError>`]. The variants are intentionally
//! coarse — a docking / screening caller usually only cares about a
//! handful of failure modes:
//!
//! 1. Did a structure / ligand / results file fail to parse
//!    ([`DockScreenError::Parse`])?
//! 2. Is the receptor itself unusable — empty, no heavy atoms, a
//!    reactive anchor atom that doesn't exist
//!    ([`DockScreenError::InvalidReceptor`])?
//! 3. Is the ligand unusable — empty, no rotatable-bond root, a
//!    covalent attachment atom that doesn't exist
//!    ([`DockScreenError::InvalidLigand`])?
//! 4. Did the caller pass nonsense arguments — a non-positive grid
//!    spacing, an empty library, an out-of-range index
//!    ([`DockScreenError::Invalid`])?
//! 5. Did an adapter ask for an external tool that isn't on `PATH`
//!    ([`DockScreenError::ToolNotAvailable`])?
//! 6. Is this a documented capability gap awaiting deeper work
//!    ([`DockScreenError::NotYetImplemented`])?
//!
//! Use [`DockScreenError::code`] for stable log / telemetry tagging
//! and [`DockScreenError::category`] to bucket failures without
//! matching every variant. The pattern mirrors `valenx-cheminf`'s
//! `CheminfError` and `valenx-biostruct`'s `BiostructError`.

use std::fmt;

/// Errors produced by `valenx-dock-screen`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DockScreenError {
    /// A structure, ligand or results file failed to parse. `format`
    /// names the expected notation (`"pdbqt"`, `"results"`, …);
    /// `detail` is a human-readable reason surfaced verbatim in the UI.
    Parse {
        /// Notation being parsed (`"pdbqt"`, `"results"`, …).
        format: &'static str,
        /// Human-readable parse-failure reason.
        detail: String,
    },

    /// A receptor is structurally unusable for docking: no atoms, no
    /// heavy atoms, a reactive anchor or flexible-sidechain selection
    /// that names an atom that does not exist. A property of the
    /// receptor, not of a parse or a caller argument.
    InvalidReceptor {
        /// Human-readable reason.
        reason: String,
    },

    /// A ligand is structurally unusable for docking: no atoms, no
    /// torsion-tree root, a covalent attachment atom that does not
    /// exist.
    InvalidLigand {
        /// Human-readable reason.
        reason: String,
    },

    /// Caller passed an argument the algorithm cannot accept: an empty
    /// library, a non-positive count or spacing, an out-of-range
    /// index, mismatched conformer / score lengths, etc. A property of
    /// the *call*.
    Invalid {
        /// Logical parameter name (e.g. `"grid_spacing"`, `"library"`).
        what: &'static str,
        /// Human-readable reason.
        reason: String,
    },

    /// A subprocess adapter was asked to run an external tool that is
    /// not installed / not on `PATH`. `tool` is the human-readable
    /// tool name; `hint` tells the user how to install it. This is the
    /// honest failure mode of every neural-network-tool adapter — the
    /// crate never reimplements AlphaFold, ProteinMPNN, DiffDock,
    /// RELION etc., it shells out to them when present.
    ToolNotAvailable {
        /// Human-readable external-tool name (e.g. `"AlphaFold"`).
        tool: &'static str,
        /// Install / PATH hint surfaced to the user.
        hint: String,
    },

    /// Feature is part of this crate's public surface but not yet
    /// implemented as a real v1. The string identifies which algorithm
    /// the caller asked for.
    NotYetImplemented {
        /// Stable feature identifier (e.g. `"macrocycle_sampling"`).
        feature: &'static str,
    },
}

/// Coarse category for routing / display.
///
/// Stable across crate versions — switch a single `match` on this
/// rather than on the error variants.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum ErrorCategory {
    /// A structure / ligand / results file failed to parse.
    Parse,
    /// User-supplied input is wrong (bad receptor, bad ligand, bad
    /// argument).
    Input,
    /// An external tool an adapter wraps is not available.
    Environment,
    /// Capability not available in v1 (documented gap).
    Capability,
}

impl DockScreenError {
    /// Stable snake-cased error code suitable for log / telemetry
    /// tagging. Format: `"dock_screen.<sub_id>"`. Codes never change
    /// across minor versions.
    pub fn code(&self) -> &'static str {
        match self {
            DockScreenError::Parse { .. } => "dock_screen.parse",
            DockScreenError::InvalidReceptor { .. } => "dock_screen.invalid_receptor",
            DockScreenError::InvalidLigand { .. } => "dock_screen.invalid_ligand",
            DockScreenError::Invalid { .. } => "dock_screen.invalid",
            DockScreenError::ToolNotAvailable { .. } => "dock_screen.tool_not_available",
            DockScreenError::NotYetImplemented { .. } => "dock_screen.not_yet_implemented",
        }
    }

    /// Coarse category string — see [`ErrorCategory`].
    pub fn category(&self) -> &'static str {
        match self {
            DockScreenError::Parse { .. } => "parse",
            DockScreenError::InvalidReceptor { .. }
            | DockScreenError::InvalidLigand { .. }
            | DockScreenError::Invalid { .. } => "input",
            DockScreenError::ToolNotAvailable { .. } => "environment",
            DockScreenError::NotYetImplemented { .. } => "capability",
        }
    }

    /// Typed category enum (for callers that want to `match` instead of
    /// comparing the [`category`](Self::category) string).
    pub fn category_enum(&self) -> ErrorCategory {
        match self {
            DockScreenError::Parse { .. } => ErrorCategory::Parse,
            DockScreenError::InvalidReceptor { .. }
            | DockScreenError::InvalidLigand { .. }
            | DockScreenError::Invalid { .. } => ErrorCategory::Input,
            DockScreenError::ToolNotAvailable { .. } => ErrorCategory::Environment,
            DockScreenError::NotYetImplemented { .. } => ErrorCategory::Capability,
        }
    }

    /// Convenience constructor for [`DockScreenError::Parse`].
    pub fn parse(format: &'static str, detail: impl Into<String>) -> Self {
        DockScreenError::Parse {
            format,
            detail: detail.into(),
        }
    }

    /// Convenience constructor for [`DockScreenError::InvalidReceptor`].
    pub fn invalid_receptor(reason: impl Into<String>) -> Self {
        DockScreenError::InvalidReceptor {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`DockScreenError::InvalidLigand`].
    pub fn invalid_ligand(reason: impl Into<String>) -> Self {
        DockScreenError::InvalidLigand {
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`DockScreenError::Invalid`].
    pub fn invalid(what: &'static str, reason: impl Into<String>) -> Self {
        DockScreenError::Invalid {
            what,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`DockScreenError::ToolNotAvailable`].
    pub fn tool_not_available(tool: &'static str, hint: impl Into<String>) -> Self {
        DockScreenError::ToolNotAvailable {
            tool,
            hint: hint.into(),
        }
    }

    /// Convenience constructor for [`DockScreenError::NotYetImplemented`].
    pub fn not_yet(feature: &'static str) -> Self {
        DockScreenError::NotYetImplemented { feature }
    }
}

impl fmt::Display for DockScreenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DockScreenError::Parse { format, detail } => {
                write!(f, "{format} parse error: {detail}")
            }
            DockScreenError::InvalidReceptor { reason } => {
                write!(f, "invalid receptor: {reason}")
            }
            DockScreenError::InvalidLigand { reason } => {
                write!(f, "invalid ligand: {reason}")
            }
            DockScreenError::Invalid { what, reason } => {
                write!(f, "invalid `{what}`: {reason}")
            }
            DockScreenError::ToolNotAvailable { tool, hint } => {
                write!(
                    f,
                    "external tool `{tool}` is not available on PATH — {hint}"
                )
            }
            DockScreenError::NotYetImplemented { feature } => {
                write!(
                    f,
                    "dock-screen feature `{feature}` is not yet implemented (v1 scaffold)"
                )
            }
        }
    }
}

impl std::error::Error for DockScreenError {}

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, DockScreenError>;

/// Bridge `valenx_dock::DockError` into [`DockScreenError`]. The dock
/// crate's parse / input failures land as [`DockScreenError::Parse`]
/// or [`DockScreenError::InvalidLigand`] depending on its category;
/// config and runtime failures land as [`DockScreenError::Invalid`].
impl From<valenx_dock::DockError> for DockScreenError {
    fn from(e: valenx_dock::DockError) -> Self {
        use valenx_dock::error::ErrorCategory as DC;
        let msg = e.to_string();
        match e.category() {
            DC::Input => DockScreenError::InvalidLigand { reason: msg },
            DC::Config => DockScreenError::Invalid {
                what: "dock_config",
                reason: msg,
            },
            DC::Runtime | DC::Internal => DockScreenError::Invalid {
                what: "dock_runtime",
                reason: msg,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn code_and_category_match_variants() {
        let err = DockScreenError::parse("pdbqt", "missing ROOT");
        assert_eq!(err.code(), "dock_screen.parse");
        assert_eq!(err.category(), "parse");
        assert_eq!(err.category_enum(), ErrorCategory::Parse);

        let err = DockScreenError::invalid_receptor("no heavy atoms");
        assert_eq!(err.code(), "dock_screen.invalid_receptor");
        assert_eq!(err.category(), "input");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = DockScreenError::invalid_ligand("no root group");
        assert_eq!(err.code(), "dock_screen.invalid_ligand");
        assert_eq!(err.category_enum(), ErrorCategory::Input);

        let err = DockScreenError::invalid("grid_spacing", "must be positive");
        assert_eq!(err.code(), "dock_screen.invalid");
        assert_eq!(err.category(), "input");

        let err = DockScreenError::tool_not_available("AlphaFold", "install ColabFold");
        assert_eq!(err.code(), "dock_screen.tool_not_available");
        assert_eq!(err.category(), "environment");
        assert_eq!(err.category_enum(), ErrorCategory::Environment);

        let err = DockScreenError::not_yet("macrocycle_sampling");
        assert_eq!(err.code(), "dock_screen.not_yet_implemented");
        assert_eq!(err.category(), "capability");
        assert_eq!(err.category_enum(), ErrorCategory::Capability);
    }

    #[test]
    fn display_is_informative() {
        let msg = DockScreenError::parse("pdbqt", "bad atom").to_string();
        assert!(msg.contains("pdbqt"), "got: {msg}");
        assert!(msg.contains("bad atom"), "got: {msg}");

        let msg =
            DockScreenError::tool_not_available("DiffDock", "pip install diffdock").to_string();
        assert!(msg.contains("DiffDock"), "got: {msg}");
        assert!(msg.contains("pip install diffdock"), "got: {msg}");
    }

    #[test]
    fn error_trait_object() {
        let err: Box<dyn std::error::Error> = Box::new(DockScreenError::invalid("x", "y"));
        assert!(err.to_string().contains('x'));
    }

    #[test]
    fn dock_error_bridges_by_category() {
        // A dock input error → InvalidLigand.
        let de: DockScreenError = valenx_dock::DockError::NoRoot.into();
        assert_eq!(de.code(), "dock_screen.invalid_ligand");
        // A dock config error → Invalid.
        let de: DockScreenError = valenx_dock::DockError::BadNumModes(0).into();
        assert_eq!(de.code(), "dock_screen.invalid");
    }
}
