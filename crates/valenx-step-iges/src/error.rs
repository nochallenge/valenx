//! Error taxonomy for the STEP / IGES bridge.
//!
//! Every public surface returns [`StepIgesError`]. Variants are
//! deliberately granular so the desktop shell can route failures to a
//! sensible toast (parse vs IO vs unsupported-feature).

use thiserror::Error;

/// All errors produced by `valenx-step-iges`.
#[derive(Debug, Error)]
pub enum StepIgesError {
    /// Caller passed a [`valenx_cad::Solid`] that had no faces / no
    /// shells — nothing to write.
    #[error("solid is empty (no shells / faces)")]
    EmptySolid,

    /// File-level parse failure: the bytes were not a well-formed STEP
    /// or IGES document.
    #[error("parse error: {0}")]
    ParseError(String),

    /// Filesystem failure (read or write).
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    /// Caller asked for a format / feature that this v1 implementation
    /// does not yet handle. The shipped surface intentionally targets
    /// the most common ~80% of CAD interchange — surface this error so
    /// the user knows what to do instead (export to STL / OBJ, file a
    /// feature request).
    #[error("unsupported in {format}: {feature}")]
    Unsupported {
        /// Format identifier — `"STEP"` or `"IGES"`.
        format: &'static str,
        /// Specific feature that wasn't handled (e.g. `"NURBS surface
        /// with rational weights"`).
        feature: String,
    },

    /// Caller tried to serialise a mesh-backed [`valenx_cad::Solid`]
    /// (the kind produced by `valenx-fillet`'s mesh-domain pipeline).
    /// STEP and IGES both want BRep topology to round-trip faces; a
    /// triangle soup is not enough. Re-export the upstream feature
    /// before the fillet, or convert to STL for triangle-mesh
    /// interchange.
    #[error("mesh-backed solid is not exportable to {format}: {reason}")]
    MeshBackedSolidNotExportable {
        /// `"STEP"` or `"IGES"`.
        format: &'static str,
        /// One-line user-facing explanation surfaced in the toast.
        reason: String,
    },

    /// Round-4 DoS hardening: a malicious IGES file advertises a list
    /// count larger than [`crate::iges::MAX_IGES_LIST_LEN`]. Without the
    /// guard, `Vec::with_capacity(usize::MAX)` would either OOM the
    /// process or trigger a panic on `Vec::with_capacity`'s allocation
    /// limit.
    #[error("iges list too large: {count} > {max} (DoS guard)")]
    ListTooLarge {
        /// Count the file advertised.
        count: usize,
        /// Cap enforced by the parser.
        max: usize,
    },

    /// Round-4 DoS hardening: arithmetic combining IGES list-count
    /// fields overflowed `usize`. Surfaces NURBS surface entries
    /// (`(k1 + m1 + 2)` for the knot vector, `(k1+1) * (k2+1)` for the
    /// control-point grid) with degenerate values designed to wrap.
    #[error("iges arithmetic overflow during list-size computation")]
    ArithmeticOverflow,

    /// Round-9 DoS hardening: the file's metadata reports a size
    /// larger than [`crate::MAX_CAD_INTERCHANGE_FILE_BYTES`]. Without
    /// the guard, `fs::read_to_string` would happily allocate and
    /// read a multi-GB file before the parser saw any of it. 256 MiB
    /// is generous for real-world STEP/IGES (production assemblies
    /// rarely exceed 50 MiB) while staying well below a hostile
    /// `cat /dev/zero > big.step` denial of service.
    #[error("{format} file too large: {size} bytes > {cap} cap (DoS guard)")]
    FileTooLarge {
        /// `"STEP"` or `"IGES"` (or `"AP242"`).
        format: &'static str,
        /// Observed size in bytes.
        size: u64,
        /// Configured cap in bytes.
        cap: u64,
    },
}

/// Coarse category for UI/LLM routing — mirrors the pattern used by
/// `valenx-feature-tree::ErrorCategory`.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ErrorCategory {
    /// Bad input data (malformed file, unsupported topology).
    Input,
    /// Transient filesystem failure (disk full, permission, etc.).
    Runtime,
    /// Caller asked for a feature this version does not implement.
    Unsupported,
}

impl StepIgesError {
    /// Stable kebab-cased identifier; new variants get new codes.
    pub fn code(&self) -> &'static str {
        match self {
            StepIgesError::EmptySolid => "step_iges.empty_solid",
            StepIgesError::ParseError(_) => "step_iges.parse_error",
            StepIgesError::Io(_) => "step_iges.io",
            StepIgesError::Unsupported { .. } => "step_iges.unsupported",
            StepIgesError::MeshBackedSolidNotExportable { .. } => {
                "step_iges.mesh_backed_solid_not_exportable"
            }
            StepIgesError::ListTooLarge { .. } => "step_iges.list_too_large",
            StepIgesError::ArithmeticOverflow => "step_iges.arithmetic_overflow",
            StepIgesError::FileTooLarge { .. } => "step_iges.file_too_large",
        }
    }

    /// Routing category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            StepIgesError::EmptySolid
            | StepIgesError::ParseError(_)
            | StepIgesError::MeshBackedSolidNotExportable { .. }
            | StepIgesError::ListTooLarge { .. }
            | StepIgesError::ArithmeticOverflow
            | StepIgesError::FileTooLarge { .. } => ErrorCategory::Input,
            StepIgesError::Io(_) => ErrorCategory::Runtime,
            StepIgesError::Unsupported { .. } => ErrorCategory::Unsupported,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_variant_has_stable_code_and_category() {
        let cases: Vec<(StepIgesError, &'static str, ErrorCategory)> = vec![
            (
                StepIgesError::EmptySolid,
                "step_iges.empty_solid",
                ErrorCategory::Input,
            ),
            (
                StepIgesError::ParseError("bad header".into()),
                "step_iges.parse_error",
                ErrorCategory::Input,
            ),
            (
                StepIgesError::Unsupported {
                    format: "STEP",
                    feature: "AP242".into(),
                },
                "step_iges.unsupported",
                ErrorCategory::Unsupported,
            ),
            (
                StepIgesError::MeshBackedSolidNotExportable {
                    format: "STEP",
                    reason: "no BRep topology".into(),
                },
                "step_iges.mesh_backed_solid_not_exportable",
                ErrorCategory::Input,
            ),
        ];
        for (err, expected_code, expected_cat) in cases {
            assert_eq!(err.code(), expected_code, "wrong code for {err:?}");
            assert_eq!(err.category(), expected_cat, "wrong category for {err:?}");
        }
    }

    #[test]
    fn io_error_wraps_std_io() {
        let inner = std::io::Error::new(std::io::ErrorKind::NotFound, "missing.step");
        let err = StepIgesError::from(inner);
        assert_eq!(err.code(), "step_iges.io");
        assert_eq!(err.category(), ErrorCategory::Runtime);
    }
}
